package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Random;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * The X-drop band of {@link GappedAligner#semiGapped} written as a uniform recurrence, which is what a row-vectorised
 * kernel has to compute (step 8 of the acceleration plan). {@link GappedAligner.Rows} follows BLAST's
 * {@code Blast_SemiGappedAlign} exactly, and there the two gap values are <em>frozen</em> at a cell the X-drop test
 * drops: {@code gapRow -= ext} and {@code gapCol -= ext} are both inside the {@code else}, and the dropped cell's
 * {@code dp[2 * bi + 1]} keeps the previous row's value. The band's left end also advances, its right end moves, and
 * the drop test's floor tightens <em>within</em> a row as the best score improves. Taken literally, every cell of the
 * row would then depend on the one before it through the branch as well as through the arithmetic.
 * <p>
 * It does not, and that is what this test pins. Two claims:
 * <ol>
 * <li><b>The freezes are unobservable.</b> At a dropped cell {@code score = max(diag, gapCol, gapRow)} is below the
 * floor, so both gap values are too; the floor never falls; and a value below the floor can only ever win a maximum
 * against values below the floor, so every score it produces is dropped in turn and no branch decision sees it. A
 * kernel may therefore decay both gaps at <em>every</em> cell.</li>
 * <li><b>The floor does not enter the arithmetic at all.</b> With the freezes gone, the row's scores no longer depend
 * on the drop decisions: the only coupling along a row is {@code gapRow}, and the diagonal comes from the row above.
 * So a kernel computes the whole row floor-free, and the X-drop is a mask applied afterwards, over a running maximum
 * that is itself a prefix maximum.</li>
 * </ol>
 * Together with {@link RowScanTest}'s two identities this leaves nothing serial in a row but the band's bookkeeping:
 * the leading run of dropped cells, the last surviving cell, and the short growth loop at the right end.
 */
class UniformBandTest {
	/**
	 * The no-gap sentinel of the uniform form. Not {@link GappedAligner#MININT}: the prefix-maximum form of the gap
	 * along the row adds {@code j * ext} to it ({@link RowScanTest}), which overflows. Far enough below any reachable
	 * score that it can never win a maximum.
	 */
	private static final int FLOOR = -(1 << 28);

	/** What one {@link GappedAligner#semiGapped} call leaves behind: its score and the best cell's offsets. */
	private record Band(int score, int aOffset, int bOffset) {
	}

	/**
	 * {@link GappedAligner#semiGapped} with the two claims applied: both gaps decay at every cell, the row's scores
	 * are computed without the floor, and the X-drop is a mask over them. Deliberately a second implementation rather
	 * than a transformation of the first, and written in the shape a vector kernel takes: a pass that reads the row
	 * above ({@code h0}), a pass that is the gap along the row, and a pass that is the bookkeeping.
	 *
	 * @param scores the score matrix, {@code Blosum62.SIZE} wide, as {@code GappedAligner.matrix.scores}
	 */
	private static Band uniform(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop,
			int gapOpen, int gapExtend, int[] scores) {
		final int n28 = Blosum62.SIZE, goe = gapOpen + gapExtend, ext = gapExtend;
		xDrop = Math.max(xDrop, goe);
		if (n <= 0 || m <= 0) {
			return new Band(0, 0, 0);
		}
		// B's letters in the order the row loop reads them, as GappedAligner.bLetters.
		final byte[] bs = new byte[n + 2];
		if (!reverse) {
			System.arraycopy(b, bBase + 1, bs, 0, Math.min(n + 1, b.length - bBase - 1));
		} else {
			for (int k = 0, end = Math.min(n + 1, bBase + n); k < end; k++) {
				bs[k] = b[bBase + n - 1 - k];
			}
		}
		// The row above: h the best score of each column, e its gap down the column (dp[2 * i] and dp[2 * i + 1]).
		final int[] h = new int[n + 2], e = new int[n + 2], h0 = new int[n + 2];
		int score = -goe;
		h[0] = 0;
		e[0] = -goe;
		int i;
		for (i = 1; i <= n; i++) {
			if (score < -xDrop) {
				break;
			}
			h[i] = score;
			e[i] = score - goe;
			score -= gapExtend;
		}
		int bSize = i, bestScore = 0, first = 0, aOffset = 0, bOffset = 0, floor = -xDrop;
		for (int ai = 1; ai <= m; ai++) {
			final int letter = reverse ? a[aBase + m - ai] : a[aBase + ai];
			// The row above, consumed: H0 = max(diagonal, gap down the column). The diagonal of the band's first cell
			// is the cell above-left, which an earlier row dropped, so there is none. Both h and e are read here and
			// nowhere else, which is why the passes below may write them in place.
			for (int j = first; j < bSize; j++) {
				final int letterB = bs[j - 1 < 0 ? 0 : j - 1] & 31;
				final int diag = j == first ? FLOOR : h[j - 1] + (letterB < n28 ? scores[letter * n28 + letterB] : 0);
				h0[j] = Math.max(diag, e[j]);
			}
			// The gap along the row, uniform: no freeze, and no floor. G(j) = max(G(j - 1) - ext, H0(j - 1) - goe),
			// which RowScanTest pins as the prefix maximum P(j - 1) - goe - (j - 1) * ext that a kernel computes.
			int gapRow = FLOOR;
			for (int j = first; j < bSize; j++) {
				final int hj = Math.max(h0[j], gapRow);
				e[j] = Math.max(hj - goe, e[j] - ext);
				gapRow = Math.max(hj - goe, gapRow - ext);
				h[j] = hj;
			}
			// The X-drop, as a mask over the scores just computed, against a floor that is a running maximum.
			int newFirst = first, last = first;
			boolean leading = true;
			for (int j = first; j < bSize; j++) {
				final int hj = h[j];
				if (hj < floor) {
					// Dropped. At the band's leading edge the cell leaves the band; inside it, its score is killed so
					// that it cannot feed a diagonal into the next row.
					if (leading) {
						newFirst = j + 1;
					} else {
						h[j] = FLOOR;
					}
				} else {
					leading = false;
					last = j;
					if (hj > bestScore) {
						bestScore = hj;
						floor = hj - xDrop;
						aOffset = ai;
						bOffset = j;
					}
				}
			}
			first = newFirst;
			if (first == bSize) {
				break;
			}
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= floor && bSize <= n) {
					h[bSize] = gapRow;
					e[bSize] = gapRow - goe;
					gapRow -= ext;
					bSize++;
				}
			}
			if (bSize <= n) {
				h[bSize] = FLOOR;
				e[bSize] = FLOOR;
				bSize++;
			}
		}
		return new Band(bestScore, aOffset, bOffset);
	}

	/** BLOSUM62 as {@code GappedAligner.matrix.scores}. */
	private static int[] blosum62() {
		final int n = Blosum62.SIZE;
		final int[] scores = new int[n * n];
		for (int x = 0; x < n; x++) {
			for (int y = 0; y < n; y++) {
				scores[x * n + y] = Blosum62.score(x, y);
			}
		}
		return scores;
	}

	@Test
	void theUniformRecurrenceScoresEveryBandAsBlastDoes() {
		final Random random = new Random(7);
		final int[] blosum62 = blosum62();
		int cells = 0, extensions = 0;
		// BLAST's costs and the ones the composition-based stage scales by 32, and cheap gaps with a tight X-drop:
		// there the band dies and revives constantly, which is what makes the floor's tightening within a row visible
		// in the result at all (at 11/1 it flips drop decisions but never the score).
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 }, { 1, 1 }, { 7, 3 },
				{ 11 * 32, 32 } }) {
			final GappedAligner aligner = new GappedAligner(costs[0], costs[1]);
			final boolean scaled = costs[1] > 1 && costs[0] > 100;
			// Half the rounds with another matrix, as the composition-based stage sets: that takes the row loop that
			// reads the matrix field rather than the copy with BLOSUM62 as constants.
			final int n = Blosum62.SIZE;
			final int[] other = new int[n * n];
			for (int x = 0; x < n; x++) {
				for (int y = 0; y < n; y++) {
					other[x * n + y] = (scaled ? 32 : 1) * Blosum62.score(x, y) + random.nextInt(5) - 2;
				}
			}
			int[] scores = scaled ? other : blosum62;
			if (scaled) {
				aligner.setMatrix(ScoreMatrix.of(other));
			}
			for (int round = 0; round < 400; round++) {
				if (round == 200 && !scaled) {
					scores = other;
					aligner.setMatrix(ScoreMatrix.of(other));
				}
				final int len = 40 + random.nextInt(260);
				final byte[] query = new byte[len + 1];
				for (int k = 1; k <= len; k++) {
					query[k] = (byte) (1 + random.nextInt(24));
				}
				// The subject: the query with substitutions, so that the band wanders and does not simply die.
				final byte[] subject = new byte[len + 2];
				for (int k = 1; k <= len; k++) {
					subject[k] = random.nextInt(4) == 0 ? (byte) (1 + random.nextInt(24)) : query[k];
				}
				// Shapes that exercise the ends: both directions, unequal lengths, and X-drops from tight (the band
				// dies and the left end chases the right) to loose (the band grows to the whole subject).
				final boolean reverse = random.nextBoolean();
				final int m = 1 + random.nextInt(len);
				final int bn = 1 + random.nextInt(len - 1);
				final int xDrop = (scaled ? 32 : 1) * (1 + random.nextInt(120));
				final int score = aligner.semiGapped(query, 0, m, subject, 0, bn, reverse, xDrop);
				final Band band = uniform(query, 0, m, subject, 0, bn, reverse, xDrop, costs[0], costs[1], scores);
				final String where = "costs " + costs[0] + "/" + costs[1] + " round " + round + " m " + m + " n " + bn
						+ (reverse ? " reverse" : " forward") + " xDrop " + xDrop;
				assertEquals(score, band.score(), where);
				assertEquals(aligner.aOffset + "," + aligner.bOffset, band.aOffset() + "," + band.bOffset(), where);
				cells += m * bn;
				extensions++;
			}
		}
		System.out.println("uniform band: " + extensions + " extensions, " + cells + " cells of band if full");
	}
}
