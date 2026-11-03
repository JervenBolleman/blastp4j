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

	private ScoreMatrix(int[] scores) {
		this.scores = scores;
	}

	/** A copy of {@code m[query][subject]} (28 x 28). */
	public static ScoreMatrix of(int[][] m) {
		final int[] s = new int[N * N];
		for (int a = 0; a < N; a++) {
			System.arraycopy(m[a], 0, s, a * N, N);
		}
		return new ScoreMatrix(s);
	}

	public int score(int query, int subject) {
		return scores[query * N + subject];
	}
}
