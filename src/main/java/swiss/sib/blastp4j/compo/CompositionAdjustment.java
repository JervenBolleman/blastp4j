package swiss.sib.blastp4j.compo;

import swiss.sib.blastp4j.stats.UngappedKarlin;

/**
 * Compositionally adjusted score matrices ({@code composition_adjustment.c} and {@code compo_mode_condition.c} of
 * BLAST+): for a query and a subject composition, decide how to adjust BLOSUM62 (mode 2, "conditional
 * compositional score matrix adjustment") and compute the adjusted, scaled integer matrix.
 * <p>
 * Matrices are in NCBIstdaa order, 28 x 28, row = query letter, column = subject letter.
 */
public final class CompositionAdjustment {
	/** Size of the NCBIstdaa alphabet. */
	public static final int ALPHSIZE = 28;
	static final int NUM_TRUE_AA = 20;
	static final int SCORE_MIN = Short.MIN_VALUE; // COMPO_SCORE_MIN = INT2_MIN

	static final int B = 2, C = 3, D = 4, E = 5, I = 9, L = 11, N = 13, Q = 15, X = 21, Z = 23, U = 24, STOP = 25,
			O = 26, J = 27;

	/** NCBIstdaa codes of the 20 true amino acids. */
	static final int[] TRUE_CHAR_POSITIONS = { 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
			22 };
	/** NCBIstdaa code to the ARNDCQEGHILKMFPSTWYV index, or -1. */
	static final int[] ALPHA_CONVERT = { -1, 0, -1, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, -1,
			18, -1, -1, -1, -1, -1 };

	private static final int LAMBDA_ITERATION_LIMIT = 100;
	private static final double LAMBDA_ERROR_TOLERANCE = 0.0000001;
	private static final double LAMBDA_FUNCTION_TOLERANCE = 0.00001;
	private static final double COMPO_ADJUST_ERR_TOLERANCE = 0.00000001;
	private static final int COMPO_ADJUST_ITERATION_LIMIT = 2000;
	/** Relative entropy of BLOSUM62. */
	private static final double FIXED_RE_BLOSUM62 = 0.44;
	private static final double MAXIMUM_X_SCORE = -1.0;
	private static final double LAMBDA_RATIO_LOWER_BOUND = 0.5;
	/** Pseudocounts for the relative-entropy based adjustment. */
	public static final int RE_PSEUDOCOUNTS = 20;

	/** How the matrix was adjusted ({@code EMatrixAdjustRule}). */
	public enum Rule {
		DONT_ADJUST, SCALE_OLD_MATRIX, UNCONSTRAINED_REL_ENTROPY, REL_ENTROPY_OLD_MATRIX_NEW_CONTEXT,
		REL_ENTROPY_OLD_MATRIX_OLD_CONTEXT, USER_SPECIFIED_REL_ENTROPY
	}

	/** Letter probabilities of a sequence ({@code Blast_AminoAcidComposition}). */
	public record Composition(double[] prob, int numTrueAminoAcids) {
	}

	/**
	 * {@code Blast_ReadAaComposition}: probabilities of the true amino acids in {@code seq[from .. from + length)};
	 * selenocysteine counts as cysteine.
	 */
	public static Composition composition(byte[] seq, int from, int length) {
		double[] prob = new double[ALPHSIZE];
		int num = 0;
		for (int i = from; i < from + length; i++) {
			int c = seq[i];
			if (ALPHA_CONVERT[c] >= 0 || c == U) {
				prob[c]++;
				num++;
			}
		}
		if (prob[U] > 0.0) {
			prob[C] += prob[U];
			prob[U] = 0.0;
		}
		if (num > 0) {
			for (int i = 0; i < ALPHSIZE; i++) {
				prob[i] /= num;
			}
		}
		return new Composition(prob, num);
	}

	/** The scaled starting matrix and its frequency ratios ({@code Blast_MatrixInfo}). */
	public static final class MatrixInfo {
		/** 29 rows of {@link #ALPHSIZE}, row-major: the last is all SCORE_MIN. */
		final int[] startMatrix = new int[(ALPHSIZE + 1) * ALPHSIZE];
		final double[][] startFreqRatios = new double[ALPHSIZE + 1][ALPHSIZE];
		final double ungappedLambda;

		/**
		 * {@code s_MatrixInfoInit} for BLOSUM62.
		 *
		 * @param idealLambda the ungapped lambda of BLOSUM62 with standard composition ({@code kbp_ideal})
		 * @param scale       the scale of the search (32)
		 */
		public MatrixInfo(double idealLambda, double scale) {
			for (int i = 0; i < ALPHSIZE; i++) {
				startMatrix[ALPHSIZE * ALPHSIZE + i] = SCORE_MIN;
				startFreqRatios[ALPHSIZE][i] = SCORE_MIN;
			}
			ungappedLambda = idealLambda / scale;
			for (int i = 0; i < ALPHSIZE; i++) {
				System.arraycopy(Blosum62Data.FREQ_RATIOS[i], 0, startFreqRatios[i], 0, ALPHSIZE);
			}
			// Blast_Int4MatrixFromFreq
			for (int i = 0; i < ALPHSIZE; i++) {
				double[] row = startFreqRatios[i].clone();
				freqRatioToScore(row, ungappedLambda);
				for (int j = 0; j < ALPHSIZE; j++) {
					startMatrix[i * ALPHSIZE + j] = round(row[j]);
				}
			}
		}

		/** The starting (unadjusted) matrix, scaled: {@link #ALPHSIZE} x {@link #ALPHSIZE}, row-major. */
		public int[] startMatrix() {
			return java.util.Arrays.copyOf(startMatrix, ALPHSIZE * ALPHSIZE);
		}
	}

	/** {@code Blast_FreqRatioToScore} on one row. */
	private static void freqRatioToScore(double[] row, double lambda) {
		for (int j = 0; j < row.length; j++) {
			row[j] = 0.0 == row[j] ? SCORE_MIN : Math.log(row[j]) / lambda;
		}
	}

	/** {@code Nint}: nearest integer, halves away from zero; values below INT_MIN become INT_MIN. */
	static int round(double x) {
		if (x < Integer.MIN_VALUE) {
			return Integer.MIN_VALUE;
		}
		x += x >= 0. ? 0.5 : -0.5;
		return (int) (long) x;
	}

	/** {@code s_GatherLetterProbs}: NCBIstdaa to the 20-letter alphabet. */
	static double[] gather(double[] stdProbs) {
		double[] out = new double[NUM_TRUE_AA];
		for (int c = 0; c < ALPHSIZE; c++) {
			if (ALPHA_CONVERT[c] != -1) {
				out[ALPHA_CONVERT[c]] = stdProbs[c];
			}
		}
		return out;
	}

	/** Workspace with BLOSUM62's joint probabilities ({@code Blast_CompositionWorkspace}). */
	public static final class Workspace {
		final double[][] matB = new double[NUM_TRUE_AA][NUM_TRUE_AA];
		final double[] firstStandardFreq = new double[NUM_TRUE_AA];
		final double[] secondStandardFreq = new double[NUM_TRUE_AA];
		final double[] matFinal = new double[NUM_TRUE_AA * NUM_TRUE_AA];

		public Workspace() {
			for (int i = 0; i < NUM_TRUE_AA; i++) {
				for (int j = 0; j < NUM_TRUE_AA; j++) {
					matB[i][j] = Blosum62Data.JOINT_PROBS[i][j];
					firstStandardFreq[i] += matB[i][j];
					secondStandardFreq[j] += matB[i][j];
				}
			}
		}
	}

	/** Result of {@link #adjustScores}. */
	public record Adjustment(int status, Rule rule, double lambdaRatio) {
	}

	/**
	 * {@code Blast_AdjustScores} with composition mode 2 and no composition p-value: fills {@code matrix} (28 x 28,
	 * row-major, scaled) with the adjusted matrix. Status 0 on success, 1 when no adjustment is possible (the matrix is then
	 * unchanged).
	 */
	public static Adjustment adjustScores(int[] matrix, Composition query, int queryLength, Composition subject,
			int subjectLength, MatrixInfo info, Workspace ws) {
		if (query.numTrueAminoAcids() == 0 || subject.numTrueAminoAcids() == 0) {
			return new Adjustment(1, Rule.DONT_ADJUST, 0);
		}
		double[] queryProbs = gather(query.prob()), matchProbs = gather(subject.prob());
		Rule rule = chooseRule(queryLength, subjectLength, queryProbs, matchProbs);
		if (rule != Rule.SCALE_OLD_MATRIX) {
			int status = compositionMatrixAdj(matrix, rule, query.numTrueAminoAcids(), subject.numTrueAminoAcids(),
					query.prob(), subject.prob(), RE_PSEUDOCOUNTS, FIXED_RE_BLOSUM62, ws, info);
			if (status <= 0) {
				return new Adjustment(status, rule, 1.0);
			}
		}
		double ratio = compositionBasedStats(matrix, info, query.prob(), subject.prob());
		return new Adjustment(0, Rule.SCALE_OLD_MATRIX, ratio);
	}

	/** {@code s_TestToApplyREAdjustmentConditional}: mode 2. */
	static Rule chooseRule(int lenQuery, int lenMatch, double[] pQuery, double[] pMatch) {
		double[] pMatrix = Blosum62Data.BACKGROUND;
		double dMMat = relativeEntropy(pMatch, pMatrix);
		double dQMat = relativeEntropy(pQuery, pMatrix);
		double dMQ = relativeEntropy(pMatch, pQuery);
		double angle = Math.acos((dMMat * dMMat + dQMat * dQMat - dMQ * dMQ) / 2.0 / dMMat / dQMat);
		angle = angle * 180 / 3.1415926543;
		double lenQ = lenQuery, lenM = lenMatch;
		double large = lenQ > lenM ? lenQ : lenM, small = lenQ > lenM ? lenM : lenQ;
		if (highPair(pQuery, lenQuery) || highPair(pMatch, lenMatch)) {
			return Rule.USER_SPECIFIED_REL_ENTROPY;
		}
		if (dMQ > 0.16 && large / small > 3.0 && angle > 70.0) {
			return Rule.SCALE_OLD_MATRIX;
		}
		return Rule.USER_SPECIFIED_REL_ENTROPY;
	}

	/** {@code s_HighPairFrequencies}. */
	private static boolean highPair(double[] p, int length) {
		if (length <= 50) {
			return false;
		}
		double max = 0, second = 0;
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			if (p[i] > second) {
				second = p[i];
				if (p[i] > max) {
					second = max;
					max = p[i];
				}
			}
		}
		return max + second > 0.4;
	}

	/** {@code Blast_GetRelativeEntropy} (the square root of the Jensen-Shannon divergence). */
	static double relativeEntropy(double[] a, double[] b) {
		double value = 0.0;
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			double temp = (a[i] + b[i]) / 2;
			if (temp > 0) {
				if (a[i] > 0) {
					value += a[i] * Math.log(a[i] / temp) / 2;
				}
				if (b[i] > 0) {
					value += b[i] * Math.log(b[i] / temp) / 2;
				}
			}
		}
		if (value < 0) {
			value = 0;
		}
		return Math.sqrt(value);
	}

	/** {@code Blast_ApplyPseudocounts}. */
	private static void applyPseudocounts(double[] probs20, int observations, double[] background, int pseudocounts) {
		double sum = 0.0;
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			sum += probs20[i];
		}
		if (sum == 0.0) {
			sum = 1.0;
		}
		double weight = (double) pseudocounts / (observations + (double) pseudocounts);
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			probs20[i] = (1.0 - weight) * probs20[i] / sum + weight * background[i];
		}
	}

	/** {@code Blast_CompositionMatrixAdj}. */
	static int compositionMatrixAdj(int[] matrix, Rule rule, int length1, int length2, double[] stdRowProbs,
			double[] stdColProbs, int pseudocounts, double specifiedRE, Workspace ws, MatrixInfo info) {
		double[] rowProbs = gather(stdRowProbs), colProbs = gather(stdColProbs);
		double desiredRe;
		switch (rule) {
		case UNCONSTRAINED_REL_ENTROPY -> desiredRe = 0.0;
		case USER_SPECIFIED_REL_ENTROPY -> desiredRe = specifiedRE;
		default -> throw new IllegalArgumentException("rule " + rule + " is not used by blastp");
		}
		applyPseudocounts(rowProbs, length1, ws.firstStandardFreq, pseudocounts);
		applyPseudocounts(colProbs, length2, ws.secondStandardFreq, pseudocounts);
		double[] q = new double[NUM_TRUE_AA * NUM_TRUE_AA];
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			System.arraycopy(ws.matB[i], 0, q, i * NUM_TRUE_AA, NUM_TRUE_AA);
		}
		int status = OptimizeTargetFreq.optimize(ws.matFinal, NUM_TRUE_AA, q, rowProbs, colProbs, desiredRe > 0.0,
				desiredRe, COMPO_ADJUST_ERR_TOLERANCE, COMPO_ADJUST_ITERATION_LIMIT);
		if (status != 0) {
			return status;
		}
		double[][] targetFreq = new double[NUM_TRUE_AA][NUM_TRUE_AA];
		for (int i = 0; i < NUM_TRUE_AA; i++) {
			System.arraycopy(ws.matFinal, i * NUM_TRUE_AA, targetFreq[i], 0, NUM_TRUE_AA);
		}
		scoresStdAlphabet(matrix, targetFreq, info.startMatrix, rowProbs, colProbs, info.ungappedLambda);
		return 0;
	}

	/** {@code s_UnpackLetterProbs} followed by {@code s_SetPairAmbigProbsToSum}. */
	private static double[] unpack(double[] probs20) {
		double[] std = new double[ALPHSIZE];
		for (int c = 0; c < ALPHSIZE; c++) {
			std[c] = ALPHA_CONVERT[c] != -1 ? probs20[ALPHA_CONVERT[c]] : 0.0;
		}
		std[B] = std[D] + std[N];
		std[Z] = std[E] + std[Q];
		std[J] = std[I] + std[L];
		return std;
	}

	/** {@code s_ScoresStdAlphabet}. */
	private static void scoresStdAlphabet(int[] matrix, double[][] targetFreq, int[] startMatrix,
			double[] rowProb, double[] colProb, double lambda) {
		double[] rowStd = unpack(rowProb), colStd = unpack(colProb);
		double[][] scores = new double[ALPHSIZE][ALPHSIZE];
		// Blast_TrueAaToStdTargetFreqs
		double sum = 0.0;
		for (int a = 0; a < NUM_TRUE_AA; a++) {
			for (int b = 0; b < NUM_TRUE_AA; b++) {
				sum += targetFreq[a][b];
			}
		}
		for (int a = 0; a < ALPHSIZE; a++) {
			if (ALPHA_CONVERT[a] < 0) {
				for (int b = 0; b < ALPHSIZE; b++) {
					scores[a][b] = 0.0;
				}
			} else {
				int ai = ALPHA_CONVERT[a];
				for (int b = 0; b < ALPHSIZE; b++) {
					scores[a][b] = ALPHA_CONVERT[b] < 0 ? 0.0 : targetFreq[ai][ALPHA_CONVERT[b]] / sum;
				}
				scores[a][B] = scores[a][D] + scores[a][N];
				scores[a][Z] = scores[a][E] + scores[a][Q];
				scores[a][J] = scores[a][I] + scores[a][L];
			}
		}
		for (int b = 0; b < ALPHSIZE; b++) {
			scores[B][b] = scores[D][b] + scores[N][b];
		}
		for (int b = 0; b < ALPHSIZE; b++) {
			scores[Z][b] = scores[E][b] + scores[Q][b];
		}
		for (int b = 0; b < ALPHSIZE; b++) {
			scores[J][b] = scores[I][b] + scores[L][b];
		}
		// Blast_CalcFreqRatios
		for (int i = 0; i < ALPHSIZE; i++) {
			if (rowStd[i] > 0) {
				for (int j = 0; j < ALPHSIZE; j++) {
					if (colStd[j] > 0) {
						scores[i][j] /= rowStd[i] * colStd[j];
					}
				}
			}
		}
		for (int i = 0; i < ALPHSIZE; i++) {
			freqRatioToScore(scores[i], lambda);
		}
		setXUOScores(scores, rowStd, colStd);
		roundInto(matrix, scores);
		stopScores(matrix, startMatrix);
	}

	private static void roundInto(int[] matrix, double[][] scores) {
		for (int p = 0; p < ALPHSIZE; p++) {
			final double[] row = scores[p];
			for (int c = 0; c < ALPHSIZE; c++) {
				matrix[p * ALPHSIZE + c] = round(row[c]);
			}
		}
	}

	/** The stop codon's row and column of {@code matrix} from the starting matrix. */
	private static void stopScores(int[] matrix, int[] startMatrix) {
		for (int i = 0; i < ALPHSIZE; i++) {
			matrix[i * ALPHSIZE + STOP] = startMatrix[i * ALPHSIZE + STOP];
		}
		System.arraycopy(startMatrix, STOP * ALPHSIZE, matrix, STOP * ALPHSIZE, ALPHSIZE);
	}

	/** {@code s_CalcAvgScore} over row {@code i} (or column, when {@code column}). */
	private static double avgScore(double[][] m, int index, boolean column, double[] probs) {
		double score = 0.0;
		for (int j = 0; j < ALPHSIZE; j++) {
			if (ALPHA_CONVERT[j] >= 0) {
				score += (column ? m[j][index] : m[index][j]) * probs[j];
			}
		}
		return score;
	}

	private static double xScore(double[][] m, int index, boolean column, double[] probs) {
		return Math.min(avgScore(m, index, column, probs), MAXIMUM_X_SCORE);
	}

	/** {@code s_SetXUOScores}. */
	private static void setXUOScores(double[][] m, double[] rowProbs, double[] colProbs) {
		double scoreXX = 0.0;
		for (int i = 0; i < ALPHSIZE; i++) {
			if (ALPHA_CONVERT[i] >= 0) {
				double avgIX = avgScore(m, i, false, colProbs);
				m[i][X] = Math.min(avgIX, MAXIMUM_X_SCORE);
				scoreXX += avgIX * rowProbs[i];
				m[X][i] = xScore(m, i, true, rowProbs);
			}
		}
		m[X][X] = Math.min(scoreXX, MAXIMUM_X_SCORE);
		m[B][X] = xScore(m, B, false, colProbs);
		m[X][B] = xScore(m, B, true, rowProbs);
		m[Z][X] = xScore(m, Z, false, colProbs);
		m[X][Z] = xScore(m, Z, true, rowProbs);
		m[J][X] = xScore(m, J, false, colProbs);
		m[X][J] = xScore(m, J, true, rowProbs);
		System.arraycopy(m[C], 0, m[U], 0, ALPHSIZE);
		for (int i = 0; i < ALPHSIZE; i++) {
			m[i][U] = m[i][C];
		}
		System.arraycopy(m[X], 0, m[O], 0, ALPHSIZE);
		for (int i = 0; i < ALPHSIZE; i++) {
			m[i][O] = m[i][X];
		}
	}

	/**
	 * {@code Blast_CompositionBasedStats}: scales the starting matrix so that its ungapped lambda in the context of the
	 * two compositions equals the standard one; returns the lambda ratio.
	 */
	static double compositionBasedStats(int[] matrix, MatrixInfo info, double[] queryProb, double[] resProb) {
		// s_GetMatrixScoreProbs
		int minScore = 0, maxScore = 0;
		for (int irow = 0; irow < ALPHSIZE; irow++) {
			for (int aa = 0; aa < NUM_TRUE_AA; aa++) {
				int v = info.startMatrix[irow * ALPHSIZE + TRUE_CHAR_POSITIONS[aa]];
				if (v < minScore && v > SCORE_MIN) {
					minScore = v;
				}
				if (v > maxScore) {
					maxScore = v;
				}
			}
		}
		double[] scoreProb = new double[maxScore - minScore + 1];
		for (int irow = 0; irow < ALPHSIZE; irow++) {
			for (int aa = 0; aa < NUM_TRUE_AA; aa++) {
				int jcol = TRUE_CHAR_POSITIONS[aa];
				int v = info.startMatrix[irow * ALPHSIZE + jcol];
				if (v >= minScore) {
					scoreProb[v - minScore] += queryProb[irow] * resProb[jcol];
				}
			}
		}
		// s_CalcLambda
		double avg = 0.0;
		for (int i = 0; i < scoreProb.length; i++) {
			avg += (minScore + i) * scoreProb[i];
		}
		double correct = UngappedKarlin.lambda(scoreProb, minScore, minScore, maxScore, avg, info.ungappedLambda);
		double ratio = correct / info.ungappedLambda;
		ratio = Math.min(1, ratio);
		ratio = Math.max(ratio, LAMBDA_RATIO_LOWER_BOUND);
		if (ratio > 0) {
			double scaledLambda = info.ungappedLambda / ratio;
			// s_ScaleSquareMatrix
			double[][] scores = new double[ALPHSIZE][];
			for (int i = 0; i < ALPHSIZE; i++) {
				scores[i] = info.startFreqRatios[i].clone();
				freqRatioToScore(scores[i], scaledLambda);
			}
			setXUOScores(scores, queryProb, resProb);
			roundInto(matrix, scores);
			stopScores(matrix, info.startMatrix);
		}
		return ratio;
	}
}
