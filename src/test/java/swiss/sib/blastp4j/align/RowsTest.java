package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * The row loops with constant gap costs: hidden copies of {@link GappedAligner.Rows} (score only) and
 * {@link GappedAligner.ExRows} (traceback), one per pair of costs.
 */
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
	void aTracebackCopyPerPairOfGapCosts() {
		// BLAST's costs, and the same scaled by 32 as the composition-based stage aligns with them.
		GappedAligner.ExRowLoops a = GappedAligner.ExRows.of(12, 1);
		GappedAligner.ExRowLoops b = GappedAligner.ExRows.of(12 * 32, 32);
		assertTrue(a.getClass().isHidden(), a.getClass().getName());
		assertTrue(b.getClass().isHidden(), b.getClass().getName());
		assertSame(a, GappedAligner.ExRows.of(12, 1));
		assertSame(b, GappedAligner.ExRows.of(12 * 32, 32));
		assertNotSame(a.getClass(), b.getClass());
	}

	@Test
	void theTracebackCopyAlignsAsTheFieldsWithScaledCosts() {
		// As CompoRedo aligns: gap costs and matrix scaled by 32, the matrix adjusted per subject.
		Random random = new Random(2);
		int[] scaled = new int[28 * 28];
		for (int x = 0; x < 28; x++) {
			for (int y = 0; y < 28; y++) {
				scaled[x * 28 + y] = 32 * swiss.sib.blastp4j.stats.Blosum62.score(x, y) + random.nextInt(17) - 8;
			}
		}
		GappedAligner copy = new GappedAligner(11 * 32, 32);
		GappedAligner.Rows plain = new GappedAligner.Rows(12 * 32, 32);
		GappedAligner fields = new GappedAligner(11 * 32, 32, plain, plain, new GappedAligner.ExRows(12 * 32, 32));
		copy.setMatrix(ScoreMatrix.of(scaled));
		fields.setMatrix(ScoreMatrix.of(scaled));
		for (int round = 0; round < 100; round++) {
			int len = 60 + random.nextInt(300);
			byte[] q = new byte[len];
			for (int i = 0; i < len; i++) {
				q[i] = (byte) (1 + random.nextInt(24));
			}
			// The subject: the query with substitutions, insertions and deletions, between sentinels.
			byte[] s = new byte[2 * len + 2];
			int n = 0;
			for (int i = 0; i < len && n < 2 * len - 2; i++) {
				int r = random.nextInt(100);
				if (r < 70) {
					s[1 + n++] = q[i];
				} else if (r < 90) {
					s[1 + n++] = (byte) (1 + random.nextInt(24));
				} else if (r < 95) {
					s[1 + n++] = q[i];
					s[1 + n++] = (byte) (1 + random.nextInt(24));
				}
			}
			byte[] subject = java.util.Arrays.copyOf(s, n + 2);
			int off = Math.min(20 + random.nextInt(len - 40), n - 1);
			copy.extendWithTraceback(q, subject, 1, n, off, off, 40 * 32);
			fields.extendWithTraceback(q, subject, 1, n, off, off, 40 * 32);
			assertEquals(fields.score, copy.score, "round " + round);
			assertEquals(fields.queryStart + "-" + fields.queryStop + "/" + fields.subjectStart + "-" + fields.subjectStop,
					copy.queryStart + "-" + copy.queryStop + "/" + copy.subjectStart + "-" + copy.subjectStop,
					"round " + round);
			assertEquals(fields.editScript.toString(), copy.editScript.toString(), "round " + round);
		}
	}

	@Test
	void theCopyAlignsAsTheFields() {
		Random random = new Random(1);
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 } }) {
			GappedAligner copy = new GappedAligner(costs[0], costs[1]);
			GappedAligner.Rows plain = new GappedAligner.Rows(costs[0] + costs[1], costs[1]);
			GappedAligner fields = new GappedAligner(costs[0], costs[1], plain, plain,
					new GappedAligner.ExRows(costs[0] + costs[1], costs[1]));
			// Half the rounds with another matrix, as the composition-based stage sets: the copy without the matrix.
			int[] other = new int[28 * 28];
			for (int x = 0; x < 28; x++) {
				for (int y = 0; y < 28; y++) {
					other[x * 28 + y] = swiss.sib.blastp4j.stats.Blosum62.score(x, y) + (x == y ? 1 : 0);
				}
			}
			for (int round = 0; round < 50; round++) {
				if (round == 25) {
					copy.setMatrix(ScoreMatrix.of(other));
					fields.setMatrix(ScoreMatrix.of(other));
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
