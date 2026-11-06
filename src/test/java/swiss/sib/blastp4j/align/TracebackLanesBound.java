package swiss.sib.blastp4j.align;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * Not a test: the upper bound of what lanes could gain in the composition-based stage's traceback. Replays recorded
 * tracebacks (query, SEG-masked subject, adjusted matrix, start, X-drop, score; recorded with a scratch hook in
 * CompoRedo) three ways: the traceback as the stage runs it, the score-only extension of the same, and the score-only
 * extensions in {@link SemiGappedLanes} with each extension's own matrix. Traceback lanes would do more than the
 * score-only lanes (a traceback byte per cell and lane, and the walk back), so unless the score-only lanes are well
 * below the scalar traceback, they cannot pay.
 * <p>
 * Usage: {@code java --add-modules jdk.incubator.vector -cp target/classes:target/test-classes
 * swiss.sib.blastp4j.align.TracebackLanesBound FILE PASSES}
 */
public class TracebackLanesBound {
	private static final int GAP_OPEN = 11 * 32, GAP_EXTEND = 32;

	public static void main(String[] args) throws IOException {
		final java.util.List<Object[]> calls = new java.util.ArrayList<>();
		try (DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(args[0]), 1 << 20))) {
			while (true) {
				final byte[] q;
				try {
					q = new byte[in.readInt()];
				} catch (EOFException e) {
					break;
				}
				in.readFully(q);
				final byte[] s = new byte[in.readInt()];
				in.readFully(s);
				final int length = in.readInt();
				final int[] m28 = new int[28 * 28], m32 = new int[1024];
				for (int x = 0; x < 28; x++) {
					for (int y = 0; y < 28; y++) {
						m28[x * 28 + y] = in.readInt();
						m32[x << 5 | y] = m28[x * 28 + y];
					}
				}
				calls.add(new Object[] { q, s, length, ScoreMatrix.of(m28), m32,
						new int[] { in.readInt(), in.readInt(), in.readInt(), in.readInt() } });
			}
		}
		final int passes = Integer.parseInt(args[1]);
		final int n = calls.size();
		int maxQ = 0;
		for (Object[] c : calls) {
			maxQ = Math.max(maxQ, ((byte[]) c[0]).length);
		}
		final GappedAligner g = new GappedAligner(GAP_OPEN, GAP_EXTEND);
		final int[] scoreOnly = new int[n];
		long tbBest = Long.MAX_VALUE, soBest = Long.MAX_VALUE, lanesBest = Long.MAX_VALUE, tlBest = Long.MAX_VALUE;
		int bad = 0;
		for (int pass = 0; pass < passes; pass++) {
			long t0 = System.nanoTime();
			for (int c = 0; c < n; c++) {
				final Object[] k = calls.get(c);
				final int[] v = (int[]) k[5];
				g.setMatrix((ScoreMatrix) k[3]);
				g.extendWithTraceback((byte[]) k[0], (byte[]) k[1], 1, (int) k[2], v[0], v[1], v[2]);
				if (pass == 0 && g.score != v[3]) {
					bad++;
				}
			}
			long t1 = System.nanoTime();
			for (int c = 0; c < n; c++) {
				final Object[] k = calls.get(c);
				final int[] v = (int[]) k[5];
				g.setMatrix((ScoreMatrix) k[3]);
				g.extendScoreOnly((byte[]) k[0], (byte[]) k[1], 1, (int) k[2], v[0], v[1], false, v[2]);
				scoreOnly[c] = g.score;
			}
			long t2 = System.nanoTime();
			final int[] lanesScore = new int[n];
			final SemiGappedLanes lanes = new SemiGappedLanes(GAP_OPEN, GAP_EXTEND, maxQ);
			for (int c = 0; c < n; c++) {
				final Object[] k = calls.get(c);
				final byte[] q = (byte[]) k[0], s = (byte[]) k[1];
				final int length = (int) k[2];
				final int[] v = (int[]) k[5], m32 = (int[]) k[4];
				final int qs = v[0], ss = v[1];
				lanes.submit(q, 0, qs + 1, s, 1, ss + 1, true, v[2], 2 * c, m32);
				if (qs + 1 < q.length && ss + 1 < length) {
					lanes.submit(q, qs, q.length - qs - 1, s, 1 + ss, length - ss - 1, false, v[2], 2 * c + 1, m32);
				}
			}
			while (lanes.step((id, score, a, b) -> lanesScore[id >> 1] += score)) {
			}
			long t3 = System.nanoTime();
			final int[] tlScore = new int[n];
			final long[] opsCount = new long[1];
			final TracebackLanes tl = new TracebackLanes(GAP_OPEN, GAP_EXTEND, maxQ);
			for (int c = 0; c < n; c++) {
				final Object[] k = calls.get(c);
				final byte[] q = (byte[]) k[0], s = (byte[]) k[1];
				final int length = (int) k[2];
				final int[] v = (int[]) k[5], m32 = (int[]) k[4];
				final int qs = v[0], ss = v[1];
				tl.submit(q, 0, qs + 1, s, 1, ss + 1, true, v[2], 2 * c, m32);
				if (qs + 1 < q.length && ss + 1 < length) {
					tl.submit(q, qs, q.length - qs - 1, s, 1 + ss, length - ss - 1, false, v[2], 2 * c + 1, m32);
				}
			}
			while (tl.step((id, score, a, b, ops) -> {
				tlScore[id >> 1] += score;
				opsCount[0] += ops.size();
			})) {
			}
			long t4 = System.nanoTime();
			if (pass == 0) {
				int tlBad = 0;
				for (int c = 0; c < n; c++) {
					tlBad += tlScore[c] != scoreOnly[c] ? 1 : 0;
				}
				System.out.printf("traceback lanes differ from score-only in %d (%d operations)%n", tlBad, opsCount[0]);
			}
			if (pass == 0) {
				int laneBad = 0;
				for (int c = 0; c < n; c++) {
					laneBad += lanesScore[c] != scoreOnly[c] ? 1 : 0;
				}
				System.out.printf("%d tracebacks, %d with another score than recorded; lanes differ from score-only in %d%n",
						n, bad, laneBad);
			}
			if (pass > 0) {
				tbBest = Math.min(tbBest, t1 - t0);
				soBest = Math.min(soBest, t2 - t1);
				lanesBest = Math.min(lanesBest, t3 - t2);
				tlBest = Math.min(tlBest, t4 - t3);
			}
		}
		System.out.printf(
				"best of %d passes: traceback %.1f ms, score-only %.1f ms, score-only lanes (%d) %.1f ms, traceback lanes %.1f ms%n",
				passes - 1, tbBest / 1e6, soBest / 1e6, SemiGappedLanes.L, lanesBest / 1e6, tlBest / 1e6);
	}
}
