package swiss.sib.blastp4j.compo;

import java.util.ArrayList;
import java.util.List;

import swiss.sib.blastp4j.align.GappedAligner;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.align.ScoreMatrix;
import swiss.sib.blastp4j.stats.Statistics;
import swiss.sib.blastp4j.stats.UngappedKarlin;

/**
 * The traceback stage with composition-based statistics ({@code Blast_RedoOneMatch} of {@code redo_alignment.c} with
 * the callbacks of {@code blast_kappa.c}): for one subject, adjust the scoring matrix to the compositions of query and
 * subject and redo the preliminary alignments with it.
 * <p>
 * All scores are scaled by 32 for precision (gap costs 352 and 32, lambda / 32). The subject is masked with SEG unless
 * the first alignment shows query and subject to be nearly identical; the query's selenocysteines are aligned as
 * cysteines.
 * <p>
 * Not thread safe: use one per thread. Instances for the same query share the immutable {@link Query}.
 */
public final class CompoRedo {
	/** Scale of scores and gap costs ({@code SCALING_FACTOR}). */
	public static final int SCALE = 32;
	private static final double LN2 = 0.69314718055994530941723212145818;
	/** Bits per aligned position above which an alignment may be nearly identical. */
	private static final double NEAR_IDENTICAL_BITS_PER_POSITION = 1.74;
	private static final int MINIMUM_LENGTH_NEAR_IDENTICAL = 50;
	/** Score difference, in bits, for an alignment to count as contained in another. */
	private static final double KAPPA_BIT_TOL = 2.0;

	/** What does not change between subjects for a query. */
	public static final class Query {
		final byte[] original;
		/** The query as aligned: selenocysteine (U) replaced by cysteine (C). */
		final byte[] data;
		final long[] words;
		final CompositionAdjustment.Composition composition;
		final CompositionAdjustment.MatrixInfo matrixInfo;

		public Query(byte[] query) {
			this.original = query;
			data = new byte[query.length];
			for (int i = 0; i < query.length; i++) {
				data[i] = query[i] != CompositionAdjustment.U ? query[i] : CompositionAdjustment.C;
			}
			words = wordArray(query);
			composition = CompositionAdjustment.composition(query, 0, query.length);
			matrixInfo = new CompositionAdjustment.MatrixInfo(UngappedKarlin.ideal().lambda(), SCALE);
		}
	}

	/**
	 * A redone subject: its alignments with an E-value within the cutoff (scores no longer scaled; possibly none), the
	 * best E-value and best scaled score, whether an alignment was computed with the matrix inherited from an earlier
	 * subject, and the last matrix adjusted for this subject (null if none).
	 */
	public record Result(List<Hsp> hsps, double bestEvalue, int bestScaledScore, boolean usedInheritedMatrix,
			ScoreMatrix finalMatrix) {
	}

	private final Query query;
	private final Statistics statistics;
	private final double lambda; // scaled
	private final int xDrop;
	private final double nearIdenticalCutoff;
	private final GappedAligner aligner;
	private final CompositionAdjustment.Workspace workspace = new CompositionAdjustment.Workspace();
	/** The adjusted matrix, row-major. */
	private final int[] matrix = new int[CompositionAdjustment.ALPHSIZE * CompositionAdjustment.ALPHSIZE];
	private final Seg seg = Seg.forCompositionBasedStatistics();

	/**
	 * @param statistics gapped statistics of the search (lambda, gap costs)
	 */
	public CompoRedo(Query query, Statistics statistics) {
		this.query = query;
		this.statistics = statistics;
		this.lambda = statistics.parameters().lambda() / SCALE;
		this.aligner = new GappedAligner(statistics.parameters().gapOpen() * SCALE,
				statistics.parameters().gapExtend() * SCALE);
		this.xDrop = (int) Math.max(25 * LN2 / lambda, (int) (25 * LN2 / statistics.parameters().lambda()));
		this.nearIdenticalCutoff = NEAR_IDENTICAL_BITS_PER_POSITION * LN2 / lambda;
	}

	/** An alignment in the composition library's terms ({@code BlastCompo_Alignment}); ends exclusive. */
	private record Align(int score, int queryStart, int queryEnd, int matchStart, int matchEnd, Hsp hsp) {
	}

	/**
	 * Redoes the preliminary HSPs of one subject ({@code residues[base .. base + length)}), best first.
	 * <p>
	 * BLAST keeps one score matrix per thread and only overwrites it when an adjustment succeeds. If the adjustment
	 * fails for the first alignment of a subject (the subject has no true amino acids left after SEG), later
	 * alignments of that subject are computed with whatever matrix the thread used last. {@code inherited} is that
	 * matrix; when null, such alignments use the matrix this instance used last, and the result reports it.
	 */
	public Result redo(List<Hsp> prelim, byte[] residues, int base, int length, double evalueCutoff,
			ScoreMatrix inherited) {
		if (inherited != null) {
			aligner.setMatrix(inherited);
		}
		final byte[] q = query.data;
		final int queryLength = query.original.length;
		List<Align> alignments = new ArrayList<>(); // head first
		byte[] subject = null;
		boolean subjectMaybeBiased = true, nearIdentical = true, oldNearIdentical = false, stale = false;
		ScoreMatrix adjustedMatrix = null; // the last adjustment for this subject
		int numAdjustments = 0;
		Align first = toAlign(prelim.get(0));
		for (int h = 0; h < prelim.size(); h++) {
			Align in = toAlign(prelim.get(h));
			if (h == 0 || subjectMaybeBiased) {
				nearIdentical = preliminaryTestNearIdentical(in, queryLength);
			}
			if (h == 0 || (subjectMaybeBiased && nearIdentical != oldNearIdentical)) {
				// s_SequenceGetProteinRange: a copy of the subject, masked with SEG when it may be biased
				subject = new byte[length + 2];
				System.arraycopy(residues, base, subject, 1, length);
				if (subjectMaybeBiased) {
					if (!nearIdentical || !testNearIdentical(subject, q, first)) {
						subjectMaybeBiased = seg.mask(subject, 1, length);
					}
				}
			}
			if (!isContained(in, alignments)) {
				boolean failed = false;
				if (h == 0 || nearIdentical != oldNearIdentical) {
					CompositionAdjustment.Composition sc = CompositionAdjustment.composition(subject, 1, length);
					CompositionAdjustment.Adjustment adj = CompositionAdjustment.adjustScores(matrix,
							query.composition, queryLength, sc, length, query.matrixInfo, workspace);
					failed = adj.status() != 0;
					if (!failed) {
						adjustedMatrix = ScoreMatrix.of(matrix);
						aligner.setMatrix(adjustedMatrix);
					}
					numAdjustments++;
				}
				if (!failed) {
					if (adjustedMatrix == null) {
						stale = true; // BLAST would use the matrix of an earlier subject here
					}
					Hsp p = in.hsp();
					aligner.extendWithTraceback(q, subject, 1, length, p.queryGappedStart, p.subjectGappedStart,
							xDrop);
					if (aligner.score >= 1) {
						Hsp n = new Hsp(aligner.score, aligner.queryStart, aligner.queryStop, aligner.subjectStart,
								aligner.subjectStop, 0, 0);
						n.editScript = aligner.editScript;
						withDistinctEnds(new Align(n.score, n.queryStart, n.queryEnd, n.subjectStart, n.subjectEnd, n),
								alignments, numAdjustments == 1);
					}
				}
			}
			oldNearIdentical = nearIdentical;
		}
		final ScoreMatrix finalMatrix = adjustedMatrix;
		if (alignments.isEmpty()) {
			return new Result(List.of(), Double.MAX_VALUE, 0, stale, finalMatrix);
		}
		// s_HSPListFromDistinctAlignments, s_HitlistReapContained, s_HitlistEvaluateAndPurge
		List<Hsp> hsps = new ArrayList<>();
		for (Align a : alignments) {
			hsps.add(a.hsp());
		}
		hsps.sort(Hsp.BY_SCORE);
		reapContained(hsps);
		double best = Double.MAX_VALUE;
		for (Hsp h : hsps) {
			h.evalue = statistics.evalueScaled(h.score, SCALE, queryLength, length);
			best = Math.min(best, h.evalue);
		}
		hsps.removeIf(h -> h.evalue > evalueCutoff);
		if (hsps.isEmpty()) {
			return new Result(List.of(), best, 0, stale, finalMatrix);
		}
		int bestScore = hsps.get(0).score;
		for (Hsp h : hsps) {
			h.score = (int) nint((double) h.score / SCALE);
		}
		return new Result(hsps, best, bestScore, stale, finalMatrix);
	}

	/** {@code BLAST_Nint}. */
	private static long nint(double x) {
		x += x >= 0. ? 0.5 : -0.5;
		return (long) x;
	}

	private static Align toAlign(Hsp p) {
		return new Align(p.score * SCALE, p.queryStart, p.queryEnd, p.subjectStart, p.subjectEnd, p);
	}

	private static boolean containedIn(int qFrom, int qTo, int q, int sFrom, int sTo, int s) {
		return qFrom <= q && qTo >= q && sFrom <= s && sTo >= s;
	}

	/** {@code s_IsContained}: inside an alignment that scores at least 2 bits more. */
	private boolean isContained(Align in, List<Align> alignments) {
		double threshold = in.score() + KAPPA_BIT_TOL * LN2 / lambda;
		for (Align a : alignments) {
			if (containedIn(a.queryStart(), a.queryEnd(), in.queryStart(), a.matchStart(), a.matchEnd(), in.matchStart())
					&& containedIn(a.queryStart(), a.queryEnd(), in.queryEnd(), a.matchStart(), a.matchEnd(), in.matchEnd())
					&& threshold <= a.score()) {
				return true;
			}
		}
		return false;
	}

	/** {@code s_WithDistinctEnds}: adds the new alignment unless one sharing an end point scores at least as well. */
	private static void withDistinctEnds(Align n, List<Align> alignments, boolean sameAdjustment) {
		for (Align a : alignments) {
			boolean matches = sameAdjustment ? sameEndPoint(n, a) : similarEndPoint(n, a);
			if (matches && n.score() <= a.score()) {
				return;
			}
		}
		alignments.removeIf(a -> sameEndPoint(n, a));
		alignments.add(0, n);
	}

	private static boolean sameEndPoint(Align n, Align a) {
		return (a.queryStart() == n.queryStart() && a.matchStart() == n.matchStart())
				|| (a.queryEnd() == n.queryEnd() && a.matchEnd() == n.matchEnd());
	}

	/** {@code s_IsSimilarEndPoint}. */
	private static boolean similarEndPoint(Align n, Align a) {
		boolean startContained = containedIn(a.queryStart(), a.queryEnd(), n.queryStart(), a.matchStart(), a.matchEnd(),
				n.matchStart());
		boolean endContained = containedIn(a.queryStart(), a.queryEnd(), n.queryEnd(), a.matchStart(), a.matchEnd(),
				n.matchEnd());
		return (startContained && n.queryStart() - n.matchStart() == a.queryStart() - a.matchStart())
				|| (endContained && n.queryEnd() - n.matchEnd() == a.queryEnd() - a.matchEnd());
	}

	/** {@code s_HitlistReapContained}: removes HSPs inside an earlier one with at least their score. */
	private static void reapContained(List<Hsp> hsps) {
		Hsp[] array = hsps.toArray(new Hsp[0]);
		for (int i = 1; i < array.length; i++) {
			Hsp h1 = array[i];
			for (int j = 0; j < i && h1 != null; j++) {
				Hsp h2 = array[j];
				if (h2 == null) {
					continue;
				}
				if (containedIn(h2.queryStart, h2.queryEnd, h1.queryStart, h2.subjectStart, h2.subjectEnd,
						h1.subjectStart)
						&& containedIn(h2.queryStart, h2.queryEnd, h1.queryEnd, h2.subjectStart, h2.subjectEnd,
								h1.subjectEnd)
						&& h1.score <= h2.score) {
					array[i] = h1 = null;
				}
			}
		}
		hsps.clear();
		for (Hsp h : array) {
			if (h != null) {
				hsps.add(h);
			}
		}
	}

	/** {@code s_preliminaryTestNearIdentical}. */
	private boolean preliminaryTestNearIdentical(Align a, int queryLength) {
		if (a.matchEnd() - a.matchStart() + 1 < Math.min(queryLength, MINIMUM_LENGTH_NEAR_IDENTICAL)) {
			return false;
		}
		int alignLength = Math.min(a.queryEnd() - a.queryStart(), a.matchEnd() - a.matchStart());
		return (double) a.score() / alignLength >= nearIdenticalCutoff;
	}

	// ---------------------------------------------------------------- near identical test (blast_kappa.c)

	/** {@code s_CreateWordArray}: hashes of the 8-letter words of the query (the last one is left 0, as in BLAST). */
	static long[] wordArray(byte[] seq) {
		final int wordSize = 8;
		if (seq.length < wordSize) {
			return null;
		}
		long mask = 0xFFFFFFFFFFL;
		long[] hashes = new long[seq.length - wordSize + 1];
		hashes[0] = hash(seq, 0, wordSize);
		for (int i = 1; i < seq.length - wordSize; i++) {
			hashes[i] = ((hashes[i - 1] << 5) & mask) + seq[i + wordSize - 1];
		}
		return hashes;
	}

	private static long hash(byte[] data, int from, int wordSize) {
		long h = 0;
		for (int k = 0; k < wordSize; k++) {
			h <<= 5;
			h += data[from + k];
		}
		return h;
	}

	/** Outputs of an extension: query and subject lengths and alignment length. */
	private static final class Ext {
		int queryLen, subjectLen, alignLen;
	}

	/** {@code s_ExtendRight} on {@code q[qo .. qo + queryLen)} and {@code s[so .. so + subjectLen)}. */
	private static int extendRight(byte[] q, int qo, int queryLen, byte[] s, int so, int subjectLen, int maxShift,
			Ext out) {
		int numIdentical = 0, qPos = 0, sPos = 0, gapsInQuery = 0, gapsInSubject = 0;
		while (qPos < queryLen && sPos < subjectLen) {
			boolean match = false;
			while (qPos < queryLen && sPos < subjectLen && q[qo + qPos] == s[so + sPos]) {
				numIdentical++;
				qPos++;
				sPos++;
			}
			for (int n = 1; n < maxShift && qPos + n + 1 < queryLen && sPos + n + 1 < subjectLen && !match; n++) {
				if (q[qo + qPos + n] == s[so + sPos + n] && q[qo + qPos + n + 1] == s[so + sPos + n + 1]) {
					qPos += n + 2;
					sPos += n + 2;
					numIdentical += 2;
					match = true;
				}
				if (!match && q[qo + qPos + n] == s[so + sPos] && q[qo + qPos + n + 1] == s[so + sPos + 1]) {
					qPos += n + 2;
					sPos += 2;
					numIdentical += 2;
					gapsInSubject += n;
					match = true;
				}
				if (!match && q[qo + qPos] == s[so + sPos + n] && q[qo + qPos + 1] == s[so + sPos + n + 1]) {
					qPos += 2;
					sPos += n + 2;
					numIdentical += 2;
					gapsInQuery += n;
					match = true;
				}
			}
			if (match) {
				continue;
			}
			break;
		}
		out.queryLen = qPos;
		out.subjectLen = sPos;
		out.alignLen = qPos > sPos ? qPos + gapsInQuery : sPos + gapsInSubject;
		return numIdentical;
	}

	/** {@code s_ExtendLeft}; adds to {@code out.alignLen}. */
	private static int extendLeft(byte[] q, int qo, int queryLen, byte[] s, int so, int subjectLen, int maxShift,
			Ext out) {
		int qPos = queryLen - 1, sPos = subjectLen - 1;
		int numIdentical = 0, gapsInQuery = 0, gapsInSubject = 0;
		while (qPos >= 0 && sPos >= 0) {
			boolean match = false;
			while (qPos > 0 && sPos > 0 && q[qo + qPos] == s[so + sPos]) {
				numIdentical++;
				qPos--;
				sPos--;
			}
			for (int n = 1; n < maxShift && qPos - n - 1 > 0 && sPos - n - 1 > 0 && !match; n++) {
				if (q[qo + qPos - n] == s[so + sPos - n] && q[qo + qPos - n - 1] == s[so + sPos - n - 1]) {
					qPos -= n + 2;
					sPos -= n + 2;
					numIdentical += 2;
					match = true;
				}
				if (!match && q[qo + qPos - n] == s[so + sPos] && q[qo + qPos - n - 1] == s[so + sPos - 1]) {
					qPos -= n + 2;
					sPos -= 2;
					numIdentical += 2;
					gapsInSubject += n;
					match = true;
				}
				if (!match && q[qo + qPos] == s[so + sPos - n] && q[qo + qPos - 1] == s[so + sPos - n - 1]) {
					qPos -= 2;
					sPos -= n + 2;
					numIdentical += 2;
					gapsInQuery += n;
					match = true;
				}
			}
			if (match) {
				continue;
			}
			break;
		}
		out.queryLen = queryLen - qPos - 1;
		out.subjectLen = subjectLen - sPos - 1;
		out.alignLen += out.queryLen > out.subjectLen ? out.queryLen + gapsInQuery : out.subjectLen + gapsInSubject;
		return numIdentical;
	}

	/** {@code s_FindNumIdentical}; {@code words} are query word hashes from index {@code wo}. */
	private static int findNumIdentical(byte[] q, int qo, long[] words, int wo, int queryLen, byte[] s, int so,
			int subjectLen, int maxShift) {
		final int wordSize = 8;
		final long mask = 0xFFFFFFFFFFL;
		if (words == null || queryLen < wordSize || subjectLen < wordSize) {
			return 0;
		}
		long hash = 0;
		int queryFrom = 0, subjectFrom = 0, numIdentical = 0;
		boolean match = false;
		Ext left = new Ext(), right = new Ext();
		for (int sPos = 0; sPos < subjectLen - wordSize; sPos++) {
			if (sPos == 0 || match) {
				hash = hash(s, so + sPos, wordSize);
			} else {
				hash = ((hash << 5) & mask) + s[so + sPos + wordSize - 1];
			}
			int qPos;
			for (qPos = queryFrom; qPos < queryLen - wordSize; qPos++) {
				if (words[wo + qPos] == hash) {
					break;
				}
			}
			if (qPos < queryLen - wordSize) {
				int queryStart = qPos, subjectStart = sPos;
				match = true;
				numIdentical += wordSize;
				left.alignLen = 0;
				numIdentical += extendLeft(q, qo + queryFrom, queryStart - queryFrom, s, so + subjectFrom,
						subjectStart - subjectFrom, maxShift, left);
				numIdentical += extendRight(q, qo + queryStart + wordSize, queryLen - queryStart - wordSize, s,
						so + subjectStart + wordSize, subjectLen - subjectStart - wordSize, maxShift, right);
				queryFrom = queryStart + wordSize + right.queryLen;
				subjectFrom = subjectStart + wordSize + right.subjectLen;
				sPos = subjectFrom - 1;
			} else {
				match = false;
			}
		}
		return numIdentical;
	}

	/**
	 * {@code s_TestNearIdentical}: whether more than 95% of the aligned part of the first alignment is identical.
	 * {@code subject} holds the subject from index 1, {@code q} the query from index 0.
	 */
	private boolean testNearIdentical(byte[] subject, byte[] q, Align align) {
		final int maxShift = 8;
		int qStart = align.queryStart, qEnd = align.queryEnd - 1;
		int sStart = align.matchStart, sEnd = align.matchEnd - 1;
		int queryLen = qEnd - qStart + 1, subjectLen = sEnd - sStart + 1;
		int alignLen = Math.min(queryLen, subjectLen);
		Ext right = new Ext(), left = new Ext();
		int numIdentical = extendRight(q, qStart, queryLen, subject, 1 + sStart, subjectLen, maxShift, right);
		if (right.queryLen >= queryLen || right.subjectLen >= subjectLen) {
			return (double) numIdentical / alignLen > 0.95;
		}
		numIdentical += extendLeft(q, qStart + right.queryLen, queryLen - right.queryLen, subject,
				1 + sStart + right.subjectLen, subjectLen - right.subjectLen, maxShift, left);
		if (left.queryLen + right.queryLen >= queryLen || left.subjectLen + right.subjectLen >= subjectLen) {
			return (double) numIdentical / alignLen > 0.95;
		}
		numIdentical += findNumIdentical(q, qStart + right.queryLen, query.words, qStart + right.queryLen,
				queryLen - left.queryLen - right.queryLen, subject, 1 + sStart + right.subjectLen,
				subjectLen - left.subjectLen - right.subjectLen, maxShift);
		return (double) numIdentical / alignLen > 0.95;
	}
}
