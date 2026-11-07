package swiss.sib.blastp4j.seed;

import swiss.sib.blastp4j.compo.Blosum62Data;
import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.blastp4j.stats.Blosum62;
import swiss.sib.blastp4j.stats.UngappedKarlin;

/**
 * BLAST's 15-letter compressed amino-acid alphabet for BLOSUM62, used by {@code blastp -task blastp-fast} to look up
 * words of 5 letters ({@code SCompressedAlphabetNew} in {@code blast_stat.c}, with 15 letters and a matrix scale of
 * 100). Each compressed letter is a group of amino acids; the groups are those of Edgar's SE-B(14) alphabet:
 * {@code ST IJV LM KR EQZ A G BD P N F Y H C W}. The other codes ({@code - X U * O}) have no compressed letter.
 * <p>
 * The score of an amino acid against a compressed letter ({@code s_BuildCompressedScoreMatrix}) is the logarithm of
 * BLOSUM62's frequency ratios with the group's amino acids, weighted by their background frequencies within the
 * group, in units of 1/100 of a bit-scaled score (100 / λ, with BLOSUM62's ungapped λ), rounded.
 */
public final class CompressedAlphabet {
	/** The groups, in order: the compressed letter of a group is its position. */
	static final String GROUPS = "ST IJV LM KR EQZ A G BD P N F Y H C W";
	/** Number of compressed letters. */
	public static final int SIZE = 15;
	/** Score scale: neighbourhood thresholds are multiplied by this ({@code kMatrixScale}). */
	static final double SCALE = 100.0;
	/** BLOSUM62's ungapped λ ({@code RPSfindUngappedLambda}: the first row of {@code blosum62_values}). */
	static final double UNGAPPED_LAMBDA = 0.3176;
	/** {@code BLAST_SCORE_MIN}. */
	private static final int SCORE_MIN = Short.MIN_VALUE;

	/** NCBIstdaa code to compressed letter; {@link #SIZE} for codes without one. */
	static final byte[] COMPRESS = new byte[Blosum62.SIZE];
	/** Scores, scaled: {@code SCORES[code * SIZE + letter]}. */
	static final int[] SCORES = new int[Blosum62.SIZE * SIZE];

	static {
		// s_BuildCompressedTranslation: the members of each group, in the string's order.
		java.util.Arrays.fill(COMPRESS, (byte) SIZE);
		final int[][] members = new int[SIZE][];
		final String[] groups = GROUPS.split(" ");
		for (int g = 0; g < SIZE; g++) {
			members[g] = new int[groups[g].length()];
			for (int i = 0; i < groups[g].length(); i++) {
				final int code = Alphabet.encode(groups[g].charAt(i));
				COMPRESS[code] = (byte) g;
				members[g][i] = code;
			}
		}
		// s_GetCompressedProbs: P(amino acid | group), from the normalized Robinson and Robinson frequencies.
		final double[] background = UngappedKarlin.robinson();
		final double[] prob = new double[Blosum62.SIZE];
		for (int g = 0; g < SIZE; g++) {
			double sum = 0.;
			for (int aa : members[g]) {
				sum += background[aa];
			}
			for (int aa : members[g]) {
				prob[aa] = background[aa] / sum;
			}
		}
		// s_BuildCompressedScoreMatrix.
		final double[] ratios = Blosum62Data.freqRatios();
		final double factor = SCALE / UNGAPPED_LAMBDA;
		final double minFreq = SCORE_MIN / factor;
		for (int q = 0; q < Blosum62.SIZE; q++) {
			for (int g = 0; g < SIZE; g++) {
				double val = 0;
				for (int aa : members[g]) {
					val += ratios[q * Blosum62.SIZE + aa] * prob[aa];
				}
				val = val < 1e-8 ? minFreq : Math.log(val);
				SCORES[q * SIZE + g] = nint(val * factor);
			}
		}
	}

	private CompressedAlphabet() {
	}

	/** {@code BLAST_Nint}: nearest integer, halves away from zero. */
	static int nint(double x) {
		return (int) (long) (x >= 0 ? x + 0.5 : x - 0.5);
	}
}
