package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * {@link TracebackLanes} gives what {@link GappedAligner#alignEx} gives, extension by extension: score, end point and
 * operations, with BLOSUM62 and with a matrix of its own per extension (as the composition-based stage's).
 */
class TracebackLanesTest {
	@Test
	void theLanesAlignAsAlignEx() {
		extend(256, false);
	}

	@Test
	void theLanesAlignAsAlignExWithAMatrixEach() {
		extend(256, true);
	}

	/** With a DP row and a letter ring of 8 slots to start with, so that they move and grow often. */
	@Test
	void theLanesAlignAsAlignExWhileTheirBuffersMoveAndGrow() {
		extend(8, true);
	}

	/** BLOSUM62 scaled by {@code scale}, each score moved by up to {@code noise} either way. */
	private static int[] matrix(Random random, int scale, int noise) {
		final int[] m = new int[28 * 28];
		for (int a = 0; a < 28; a++) {
			for (int b = 0; b < 28; b++) {
				m[a * 28 + b] = swiss.sib.blastp4j.stats.Blosum62.score(a, b) * scale
						+ (noise > 0 ? random.nextInt(2 * noise + 1) - noise : 0);
			}
		}
		return m;
	}

	private static void extend(int slots, boolean matrices) {
		Random random = new Random(1);
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 }, { 11 * 32, 32 } }) {
			final int scale = costs[1] == 32 ? 32 : 1;
			if (scale == 1 && matrices && costs[0] != 11) {
				continue;
			}
			GappedAligner aligner = new GappedAligner(costs[0], costs[1]);
			TracebackLanes lanes = new TracebackLanes(costs[0], costs[1], 400, slots);
			// Related pairs (a mutated copy, with insertions and deletions, so that bands drift apart, some with a long
			// insertion, so that a band drifts far from the others) and unrelated ones, extended both ways from random
			// points, more of them than there are lanes.
			List<Object[]> expected = new ArrayList<>();
			Object[][] got = new Object[300][];
			TracebackLanes.Done done = (id, sc, a, b, ops) -> got[id] = new Object[] { sc, a, b, ops.toString() };
			for (int c = 0; c < got.length; c++) {
				int len = 20 + random.nextInt(380);
				byte[] q = new byte[len];
				for (int i = 0; i < len; i++) {
					q[i] = (byte) (1 + random.nextInt(24));
				}
				byte[] s = new byte[len + 602];
				int sl = 0;
				boolean related = random.nextInt(4) != 0;
				int block = related && random.nextInt(5) == 0 ? random.nextInt(len / 2) : -1;
				for (int i = 0; i < len && sl < len + 600; i++) {
					if (i == block) {
						for (int k = 100 + random.nextInt(300); k > 0; k--) {
							s[1 + sl++] = (byte) (1 + random.nextInt(24));
						}
					}
					int r = random.nextInt(20);
					if (related && r == 0) {
						continue; // deletion
					}
					if (related && r == 1) {
						s[1 + sl++] = (byte) (1 + random.nextInt(24)); // insertion
					}
					s[1 + sl++] = related && r > 4 ? q[i] : (byte) (1 + random.nextInt(24));
				}
				int qOff = random.nextInt(len), sOff = random.nextInt(sl);
				boolean reverse = random.nextBoolean();
				int xDrop = (20 + random.nextInt(60)) * scale;
				int aBase = reverse ? 0 : qOff, m = reverse ? qOff + 1 : len - qOff - 1;
				int bBase = reverse ? 1 : 1 + sOff, n = reverse ? sOff + 1 : sl - sOff - 1;
				if (m <= 0 || n <= 0) {
					continue;
				}
				final int[] m28 = matrix(random, scale, matrices ? 3 * scale : 0), m32 = new int[1024];
				for (int a = 0; a < 28; a++) {
					for (int b = 0; b < 28; b++) {
						m32[a << 5 | b] = m28[a * 28 + b];
					}
				}
				aligner.setMatrix(ScoreMatrix.of(m28));
				EditScript ops = new EditScript();
				int score = aligner.alignEx(q, aBase, m, s, bBase, n, reverse, xDrop, ops);
				expected.add(new Object[] { c, score, aligner.aOffset, aligner.bOffset, ops.toString() });
				lanes.submit(q, aBase, m, s, bBase, n, reverse, xDrop, c, m32);
				if (random.nextInt(3) == 0) {
					lanes.step(done);
				}
			}
			while (lanes.step(done)) {
			}
			for (Object[] e : expected) {
				final int c = (int) e[0];
				assertEquals(e[1], got[c][0], "score of " + c);
				assertEquals(e[2], got[c][1], "aOffset of " + c);
				assertEquals(e[3], got[c][2], "bOffset of " + c);
				assertEquals(e[4], got[c][3], "operations of " + c);
			}
		}
	}
}
