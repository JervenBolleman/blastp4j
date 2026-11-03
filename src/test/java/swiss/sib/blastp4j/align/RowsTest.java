package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

/** The row loops with constant gap costs: a hidden copy of {@link GappedAligner.Rows}, one per pair of costs. */
class RowsTest {
	@Test
	void aCopyPerPairOfGapCostsAndMatrix() {
		GappedAligner.RowLoops a = GappedAligner.Rows.of(12, 1, null);
		GappedAligner.RowLoops b = GappedAligner.Rows.of(12, 1, GappedAligner.Rows.BLOSUM62);
		assertTrue(a.getClass().isHidden(), a.getClass().getName());
		assertTrue(b.getClass().isHidden(), b.getClass().getName());
		assertSame(a, GappedAligner.Rows.of(12, 1, null));
		assertSame(b, GappedAligner.Rows.of(12, 1, GappedAligner.Rows.BLOSUM62));
		assertNotSame(a.getClass(), b.getClass());
		assertNotSame(a.getClass(), GappedAligner.Rows.of(11, 2, null).getClass());
	}

	@Test
	void theCopyAlignsAsTheFields() {
		Random random = new Random(1);
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 } }) {
			GappedAligner copy = new GappedAligner(costs[0], costs[1]);
			GappedAligner.Rows plain = new GappedAligner.Rows(costs[0] + costs[1], costs[1]);
			GappedAligner fields = new GappedAligner(costs[0], costs[1], plain, plain);
			// Half the rounds with another matrix, as the composition-based stage sets: the copy without the matrix.
			int[][] other = new int[28][28];
			for (int x = 0; x < 28; x++) {
				for (int y = 0; y < 28; y++) {
					other[x][y] = swiss.sib.blastp4j.stats.Blosum62.score(x, y) + (x == y ? 1 : 0);
				}
			}
			for (int round = 0; round < 50; round++) {
				if (round == 25) {
					copy.setMatrix(other);
					fields.setMatrix(other);
				}
				// A random query, and a mutated copy of it as the subject, between sentinels as in the database.
				int len = 60 + random.nextInt(200);
				byte[] q = new byte[len];
				for (int i = 0; i < len; i++) {
					q[i] = (byte) (1 + random.nextInt(24));
				}
				byte[] s = new byte[len + 2];
				for (int i = 0; i < len; i++) {
					s[i + 1] = random.nextInt(5) == 0 ? (byte) (1 + random.nextInt(24)) : q[i];
				}
				int off = 20 + random.nextInt(len - 40);
				copy.extendWithTraceback(q, s, 1, len, off, off, 40);
				fields.extendWithTraceback(q, s, 1, len, off, off, 40);
				assertEquals(fields.score, copy.score);
				assertEquals(fields.editScript.toString(), copy.editScript.toString());
				copy.extendScoreOnly(q, s, 1, len, off, off, false, 40);
				fields.extendScoreOnly(q, s, 1, len, off, off, false, 40);
				assertEquals(fields.score, copy.score);
				assertEquals(fields.queryStart, copy.queryStart);
				assertEquals(fields.subjectStop, copy.subjectStop);
			}
		}
	}
}
