package swiss.sib.blastp4j.align;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * A score matrix for the gapped aligner, {@code score(query letter, subject letter)} over NCBIstdaa codes (28 x 28).
 * Immutable, so one instance can be shared between aligners and threads.
 */
public final class ScoreMatrix {
	static final int N = Blosum62.SIZE;

	/** BLOSUM62, unscaled. */
	public static final ScoreMatrix BLOSUM62;
	static {
		final int[] s = new int[N * N];
		for (int a = 0; a < N; a++) {
			for (int b = 0; b < N; b++) {
				s[a * N + b] = Blosum62.score(a, b);
			}
		}
		BLOSUM62 = new ScoreMatrix(s);
	}

	/** The scores, {@code scores[query * N + subject]}. */
	final int[] scores;
	/**
	 * The same scores 32 wide, {@code wide[query << 5 | subject]}, zero outside 28 x 28: 1024 entries, so that an index
	 * masked with {@code wide.length - 1} needs no bounds check in C2's loops. Shorts, so that the table is 2 KB rather
	 * than 4 KB in the L1D {@link GappedAligner.ExRows}' cell loops share with their band arrays and trace: every score
	 * fits, since BLAST's own floor for a matrix entry is {@code COMPO_SCORE_MIN = INT2_MIN} and the scaled scores stay
	 * far below INT2_MAX (the widest seen over a full search of the human proteome is [-32768, 533]). The loops widen
	 * each score to an int as they load it, which is free on both ISAs, so the arithmetic is unchanged.
	 */
	final short[] wide;

	private ScoreMatrix(int[] scores) {
		this.scores = scores;
		this.wide = new short[32 * 32];
		for (int a = 0; a < N; a++) {
			for (int b = 0; b < N; b++) {
				final int v = scores[a * N + b];
				if (v < Short.MIN_VALUE || v > Short.MAX_VALUE) {
					throw new IllegalArgumentException("score out of 16-bit range: " + v);
				}
				wide[(a << 5) + b] = (short) v;
			}
		}
	}

	/** A copy of {@code m[query * 28 + subject]} (28 x 28, row-major). */
	public static ScoreMatrix of(int[] m) {
		return new ScoreMatrix(java.util.Arrays.copyOf(m, N * N));
	}

	public int score(int query, int subject) {
		return scores[query * N + subject];
	}
}
