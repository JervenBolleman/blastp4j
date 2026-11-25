package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.util.Random;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * {@link RowLanes}, the row loop that computes several cells of a row at once, against
 * {@link GappedAligner.Rows}, which follows {@code Blast_SemiGappedAlign} one cell at a time. The score and the best
 * cell have to be the same in every extension: the lanes are an identity, not an approximation.
 */
class RowLanesTest {
	/** BLOSUM62 as the aligner's matrix. */
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
	void onlyForAMatrixThatFitsAByte() {
		final int[] scores = blosum62();
		org.junit.jupiter.api.Assertions.assertTrue(RowLanes.fits(scores));
		final int[] scaled = scores.clone();
		for (int k = 0; k < scaled.length; k++) {
			scaled[k] *= 32;
		}
		org.junit.jupiter.api.Assertions.assertFalse(RowLanes.fits(scaled));
		assertEquals(null, RowLanes.of(12, 1, scaled, cells(12, 1)));
		assertNotNull(RowLanes.of(12, 1, scores, cells(12, 1)));
	}

	@Test
	void theLanesScoreEveryBandAsTheCellsDo() {
		final Random random = new Random(13);
		final int[] scores = blosum62();
		int extensions = 0;
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 }, { 1, 1 }, { 7, 3 } }) {
			final GappedAligner.RowLoops lanes = RowLanes.of(costs[0] + costs[1], costs[1], scores,
					cells(costs[0] + costs[1], costs[1]));
			assertNotNull(lanes);
			final GappedAligner cells = new GappedAligner(costs[0], costs[1]);
			final GappedAligner vector = new GappedAligner(costs[0], costs[1], lanes, lanes);
			for (int round = 0; round < 400; round++) {
				final int len = 40 + random.nextInt(260);
				final byte[] query = new byte[len + 1];
				for (int k = 1; k <= len; k++) {
					query[k] = (byte) (1 + random.nextInt(24));
				}
				final byte[] subject = new byte[len + 2];
				for (int k = 1; k <= len; k++) {
					subject[k] = random.nextInt(4) == 0 ? (byte) (1 + random.nextInt(24)) : query[k];
				}
				final boolean reverse = random.nextBoolean();
				final int m = 1 + random.nextInt(len), bn = 1 + random.nextInt(len - 1);
				final int xDrop = 1 + random.nextInt(120);
				final int one = cells.semiGapped(query, 0, m, subject, 0, bn, reverse, xDrop);
				final int many = vector.semiGapped(query, 0, m, subject, 0, bn, reverse, xDrop);
				final String where = "costs " + costs[0] + "/" + costs[1] + " round " + round + " m " + m + " n " + bn
						+ (reverse ? " reverse" : " forward") + " xDrop " + xDrop;
				assertEquals(one, many, where);
				assertEquals(cells.aOffset + "," + cells.bOffset, vector.aOffset + "," + vector.bOffset, where);
				extensions++;
			}
		}
		System.out.println("row lanes: " + extensions + " extensions scored as the cells do");
	}

	/**
	 * Many small bands, where the groups of lanes are few and the band's ends move every row: the shapes that found
	 * the one real defect of the first kernel (a group's store overwrites the score the next group's diagonal reads).
	 */
	@Test
	void theLanesScoreSmallBandsAsTheCellsDo() {
		final Random random = new Random(5);
		final int[] scores = blosum62();
		final GappedAligner.RowLoops lanes = RowLanes.of(12, 1, scores, cells(12, 1));
		final GappedAligner cells = new GappedAligner(11, 1);
		final GappedAligner vector = new GappedAligner(11, 1, lanes, lanes);
		for (int round = 0; round < 200000; round++) {
			final int len = 4 + random.nextInt(26);
			final byte[] query = new byte[len + 1];
			for (int k = 1; k <= len; k++) {
				query[k] = (byte) (1 + random.nextInt(24));
			}
			final byte[] subject = new byte[len + 2];
			for (int k = 1; k <= len; k++) {
				subject[k] = random.nextInt(4) == 0 ? (byte) (1 + random.nextInt(24)) : query[k];
			}
			final boolean reverse = random.nextBoolean();
			final int m = 1 + random.nextInt(len), bn = 1 + random.nextInt(len - 1);
			final int xDrop = 1 + random.nextInt(40);
			final int one = cells.semiGapped(query, 0, m, subject, 0, bn, reverse, xDrop);
			final int many = vector.semiGapped(query, 0, m, subject, 0, bn, reverse, xDrop);
			final String where = "round " + round + " m " + m + " n " + bn + (reverse ? " reverse" : " forward")
					+ " xDrop " + xDrop;
			assertEquals(one, many, where);
			assertEquals(cells.aOffset + "," + cells.bOffset, vector.aOffset + "," + vector.bOffset, where);
		}
	}

	/** The cell-at-a-time row loop the lanes hand an extension that 16 bits cannot hold. */
	private static GappedAligner.RowLoops cells(int goe, int ext) {
		return GappedAligner.Rows.of(goe, ext, GappedAligner.Rows.BLOSUM62);
	}

	/**
	 * The two bounds that send an extension to the cells, each one taken past its limit: a band long enough for a score
	 * to pass 16 bits, and an X-drop deep enough that the sentinel would no longer behave as BLAST's MININT. The
	 * answer has to be the cells' answer either way, which is the whole point of handing it over.
	 */
	@Test
	void anExtensionTooLongOrTooDeepForSixteenBitsGoesToTheCells() {
		final Random random = new Random(29);
		final GappedAligner.RowLoops lanes = RowLanes.of(12, 1, blosum62(), cells(12, 1));
		assertNotNull(lanes);
		final GappedAligner cells = new GappedAligner(11, 1);
		final GappedAligner vector = new GappedAligner(11, 1, lanes, lanes);
		// 4,000 letters of a sequence against itself scores about 4 * 4000, past what a short holds once the prefix
		// maximum's lift is added, and min(m, n) is past the band limit.
		final int len = 4000;
		final byte[] query = new byte[len + 1];
		for (int k = 1; k <= len; k++) {
			query[k] = (byte) (1 + random.nextInt(24));
		}
		final byte[] subject = new byte[len + 2];
		System.arraycopy(query, 1, subject, 1, len);
		for (int[] shape : new int[][] { { len, len, 300 }, { 40, 40, 9000 }, { 40, 40, 100 } }) {
			final int one = cells.semiGapped(query, 0, shape[0], subject, 0, shape[1], false, shape[2]);
			final int many = vector.semiGapped(query, 0, shape[0], subject, 0, shape[1], false, shape[2]);
			final String where = "m " + shape[0] + " n " + shape[1] + " xDrop " + shape[2];
			assertEquals(one, many, where);
			assertEquals(cells.aOffset + "," + cells.bOffset, vector.aOffset + "," + vector.bOffset, where);
			org.junit.jupiter.api.Assertions.assertTrue(one > 0, where);
		}
	}
}
