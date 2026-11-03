package swiss.sib.blastp4j.align;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

/** {@link SemiGappedLanes} gives what {@link GappedAligner#semiGapped} gives, extension by extension. */
class SemiGappedLanesTest {
	@Test
	void theLanesExtendAsSemiGapped() {
		extend(256);
	}

	/** With a DP row and a letter ring of 8 slots to start with, so that they move and grow often. */
	@Test
	void theLanesExtendAsSemiGappedWhileTheirBuffersMoveAndGrow() {
		extend(8);
	}

	/** Extensions too long for the row loop's shorts, among ones that fit, are reported as semiGapped computes them. */
	@Test
	void extensionsTooLongForShortsAsSemiGapped() {
		Random random = new Random(2);
		GappedAligner aligner = new GappedAligner(11, 1);
		SemiGappedLanes lanes = new SemiGappedLanes(11, 1, 2000, 64);
		List<int[]> expected = new ArrayList<>();
		int[][] got = new int[40][];
		for (int c = 0; c < got.length; c++) {
			int m = c % 3 == 0 ? 1500 + random.nextInt(400) : 20 + random.nextInt(200);
			byte[] q = new byte[m + 1];
			byte[] s = new byte[m + 2];
			for (int i = 0; i <= m; i++) {
				q[i] = (byte) (1 + random.nextInt(24));
				s[1 + i] = random.nextInt(10) == 0 ? (byte) (1 + random.nextInt(24)) : q[i];
			}
			assertEquals(c % 3 != 0, SemiGappedLanes.fits(m, m));
			int score = aligner.semiGapped(q, 0, m, s, 1, m, false, 40);
			expected.add(new int[] { c, score, aligner.aOffset, aligner.bOffset });
			lanes.submit(q, 0, m, s, 1, m, false, 40, c);
		}
		while (lanes.step((id, sc, a, b) -> got[id] = new int[] { sc, a, b })) {
		}
		for (int[] e : expected) {
			assertEquals(e[1], got[e[0]][0], "score of " + e[0]);
			assertEquals(e[2], got[e[0]][1], "aOffset of " + e[0]);
			assertEquals(e[3], got[e[0]][2], "bOffset of " + e[0]);
		}
	}

	private static void extend(int slots) {
		Random random = new Random(1);
		for (int[] costs : new int[][] { { 11, 1 }, { 10, 1 }, { 9, 2 }, { 13, 3 } }) {
			GappedAligner aligner = new GappedAligner(costs[0], costs[1]);
			SemiGappedLanes lanes = new SemiGappedLanes(costs[0], costs[1], 400, slots);
			// Related pairs (a mutated copy, with insertions and deletions, so that bands drift apart, some with a long
			// insertion, so that a band drifts far from the others) and unrelated ones, extended both ways from random
			// points, more of them than there are lanes.
			List<int[]> expected = new ArrayList<>();
			int[][] got = new int[300][];
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
				int xDrop = 20 + random.nextInt(60);
				int aBase = reverse ? 0 : qOff, m = reverse ? qOff + 1 : len - qOff - 1;
				int bBase = reverse ? 1 : 1 + sOff, n = reverse ? sOff + 1 : sl - sOff - 1;
				if (m <= 0 || n <= 0) {
					continue;
				}
				int score = aligner.semiGapped(q, aBase, m, s, bBase, n, reverse, xDrop);
				expected.add(new int[] { c, score, aligner.aOffset, aligner.bOffset });
				lanes.submit(q, aBase, m, s, bBase, n, reverse, xDrop, c);
				if (random.nextInt(3) == 0) {
					lanes.step((id, sc, a, b) -> got[id] = new int[] { sc, a, b });
				}
			}
			while (lanes.step((id, sc, a, b) -> got[id] = new int[] { sc, a, b })) {
			}
			for (int[] e : expected) {
				assertEquals(e[1], got[e[0]][0], "score of " + e[0]);
				assertEquals(e[2], got[e[0]][1], "aOffset of " + e[0]);
				assertEquals(e[3], got[e[0]][2], "bOffset of " + e[0]);
			}
		}
	}
}
