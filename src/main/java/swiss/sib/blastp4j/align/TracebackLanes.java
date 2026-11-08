package swiss.sib.blastp4j.align;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShape;
import jdk.incubator.vector.VectorSpecies;

/**
 * Many independent {@link GappedAligner#alignEx} extensions (a score matrix each) at once, one per lane: the DP of
 * {@link SemiGappedLanes} (alignEx's cells are semiGapped's), with each cell's traceback byte as alignEx computes it,
 * and alignEx's walk back from the best cell when a lane's extension ends. Each lane gives exactly alignEx's score,
 * end point and operations.
 * <p>
 * The traceback bytes of a row step are kept per step, slot by slot with a byte per lane, from the lowest slot of the
 * bands: {@code trace[t & (steps - 1)]}, from slot {@code traceLo[t & (steps - 1)]}. An extension's rows are at most
 * {@code maxA} steps, so with more steps than that in the ring, a step's bytes are kept until every extension that
 * has a row in it has ended. Not thread safe: use one per thread.
 */
public final class TracebackLanes implements TracebackKernel {
	/**
	 * Receives the result of an extension: alignEx's score, its {@code aOffset} and {@code bOffset}, and its
	 * operations, from the best cell back to the start (valid during the call only).
	 */
	public interface Done {
		void done(int id, int score, int aOffset, int bOffset, EditScript ops);
	}

	/**
	 * The machine's preferred vectors, as SemiGappedLanes' (bin/jvm.options chooses the width), where they have the 8
	 * ints or more that the lanes run with; 256 bits elsewhere (AArch64's 128), where only the tests use the class: the
	 * traceback bytes need a byte vector of 8 x lanes bits, at least 64.
	 */
	private static final VectorSpecies<Integer> S = IntVector.SPECIES_PREFERRED.length() >= 8 ? IntVector.SPECIES_PREFERRED
			: IntVector.SPECIES_256;
	/** Lanes. */
	public static final int L = S.length();
	private static final VectorSpecies<Byte> B = VectorSpecies.of(byte.class, VectorShape.forBitSize(8 * L));
	private static final int SUB = EditScript.SUB, INS = EditScript.INS, DEL = EditScript.DEL;
	private static final int EXTEND_GAP_A = GappedAligner.EXTEND_GAP_A, EXTEND_GAP_B = GappedAligner.EXTEND_GAP_B;
	private static final int MININT = GappedAligner.MININT;
	private static final int[] MAT = GappedAligner.Rows.BLOSUM62;

	private final int goe, ext;
	/** The longest A an extension may have. */
	private final int maxA;
	/** The slot of cell (0, 0): diagonal coordinates are {@code bi - ai + off}. */
	private final int off;
	/** The DP row: slots {@code [dpBase, dpBase + dpSlots)}, slot k of lane l at {@code (k - dpBase) * L + l}. */
	private int[] dpBest, dpGap;
	private int dpBase, dpSlots;
	/**
	 * B's letters by lane, in a ring of {@code slots} positions, twice, so that a run of up to {@code slots} positions
	 * from any start is contiguous: lane l's letter for slot k in row step t is at position {@code k + t}. An extension
	 * whose row ai comes at step {@code ai + shift} has letter bi (as semiGapped reads it) at position
	 * {@code bi + shift + off + 1}.
	 */
	private int[] ring;
	/**
	 * The positions in the letter ring: more than the span of the bands of all lanes plus the cells a row's gap may
	 * add after a band and the letters written ahead.
	 */
	private int slots;
	/** The most cells a row's gap may add after a band: X-drop over the gap extension cost. */
	private int reach;
	/** A's letters as matrix rows ({@code letter << 5}) by lane: row ai at position {@code ai + shift}. */
	private final int[] aRing;
	private final int aSize;
	/** The row step. */
	private int t;

	// Per lane: the extension (as semiGapped's arguments), and semiGapped's state.
	private final int[] id = new int[L], ai = new int[L], m = new int[L], n = new int[L], first = new int[L],
			bSize = new int[L], best = new int[L], aOffset = new int[L], bOffset = new int[L], xDrop = new int[L],
			gapRow = new int[L], last = new int[L], kStart = new int[L], width = new int[L], row = new int[L],
			bBase = new int[L], shift = new int[L], written = new int[L], busy = new int[L];
	private final boolean[] reverse = new boolean[L];
	private final byte[][] subject = new byte[L][];
	private final int[] tmp = new int[L], idx = new int[L];
	/** The traceback bytes by row step, and the slot each step's start; {@code steps} entries (a power of two). */
	private final byte[][] trace;
	private final int[] traceLo;
	private final int steps;
	/** The operations of the extension being reported. */
	private final EditScript ops = new EditScript();
	/**
	 * Each lane's score matrix, 32 x 32 ({@code [a << 5 | b]}), at {@code l << 10}: BLOSUM62, or the matrix its
	 * extension was submitted with; {@code laneMatrix} is which.
	 */
	private final int[] mats = new int[L << 10];
	private final long[] laneMatrix = new long[L];
	private static final IntVector LANE_BASE = IntVector.zero(S).addIndex(1).lanewise(VectorOperators.LSHL, 10);

	// The queue of submitted extensions.
	private byte[][] qa = new byte[16][], qb = new byte[16][];
	private int[] qaBase = new int[16], qm = new int[16], qbBase = new int[16], qn = new int[16], qx = new int[16],
			qid = new int[16];
	private boolean[] qrev = new boolean[16];
	private int[][] qmat = new int[16][];
	private long[] qver = new long[16];
	private int head, size;

	/** Lanes for the given gap costs, and for extensions with A of at most {@code maxA} letters. */
	public TracebackLanes(int gapOpen, int gapExtend, int maxA) {
		this(gapOpen, gapExtend, maxA, 256);
	}

	/** With DP row and letter ring of {@code slots} slots to start with (the tests start small, to move and grow). */
	TracebackLanes(int gapOpen, int gapExtend, int maxA, int slots) {
		this.slots = slots;
		this.dpSlots = slots;
		goe = gapOpen + gapExtend;
		ext = gapExtend;
		this.maxA = maxA;
		off = maxA + 2;
		dpBest = new int[dpSlots * L];
		dpGap = new int[dpSlots * L];
		dpBase = off - dpSlots / 2;
		ring = new int[2 * slots * L];
		aSize = Integer.highestOneBit(maxA + 2) * 2;
		aRing = new int[aSize * L];
		steps = Integer.highestOneBit(maxA + 2) * 2;
		trace = new byte[steps][];
		java.util.Arrays.fill(laneMatrix, Long.MIN_VALUE);
		traceLo = new int[steps];
		for (int l = 0; l < L; l++) {
			idle(l);
		}
	}

	/** Whether these lanes take extensions with these gap costs and A of up to {@code maxA} letters. */
	public boolean fits(int gapOpen, int gapExtend, int maxA) {
		return goe == gapOpen + gapExtend && ext == gapExtend && maxA <= this.maxA;
	}

	/** The longest A an extension may have. */
	public int maxA() {
		return maxA;
	}

	/** Extensions submitted and not yet in a lane. */
	public int queued() {
		return size;
	}

	/**
	 * Submits {@code semiGapped(a, aBase, m, b, bBase, n, reverse, xDrop)}, with {@code 0 < m <= maxA} and
	 * {@code n > 0}; {@link #step} reports its result under {@code id}.
	 */
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id) {
		submit(a, aBase, m, b, bBase, n, reverse, xDrop, id, MAT, 0);
	}

	/** As {@link #submit}, with the score matrix {@code matrix[a << 5 | b]} (32 x 32) instead of BLOSUM62. */
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id,
			int[] matrix) {
		submit(a, aBase, m, b, bBase, n, reverse, xDrop, id, matrix, -1 - System.identityHashCode(matrix));
	}

	/**
	 * As {@link #submit}, with the matrix's {@code version}: a lane copies the matrix only where its version is not
	 * the one the lane has (0 is BLOSUM62), so the caller may refill one array with another matrix under a new
	 * version, once no extension it submitted with the old one is queued.
	 */
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id,
			int[] matrix, long version) {
		if (m <= 0 || n <= 0 || m > maxA) {
			throw new IllegalArgumentException("m " + m + ", n " + n);
		}
		reach = Math.max(reach, Math.max(xDrop, goe) / ext + 2);
		if (size == qa.length) {
			grow();
		}
		final int i = (head + size) % qa.length;
		qa[i] = a;
		qaBase[i] = aBase;
		qm[i] = m;
		qb[i] = b;
		qbBase[i] = bBase;
		qn[i] = n;
		qrev[i] = reverse;
		qx[i] = xDrop;
		qid[i] = id;
		qmat[i] = matrix;
		qver[i] = version;
		size++;
	}

	private void grow() {
		final int c = qa.length, c2 = 2 * c;
		final byte[][] a = new byte[c2][], b = new byte[c2][];
		final int[] aBase = new int[c2], mm = new int[c2], bBase = new int[c2], nn = new int[c2], x = new int[c2],
				ids = new int[c2];
		final boolean[] rev = new boolean[c2];
		final int[][] mat = new int[c2][];
		final long[] ver = new long[c2];
		for (int k = 0; k < size; k++) {
			final int i = (head + k) % c;
			a[k] = qa[i];
			aBase[k] = qaBase[i];
			mm[k] = qm[i];
			b[k] = qb[i];
			bBase[k] = qbBase[i];
			nn[k] = qn[i];
			rev[k] = qrev[i];
			x[k] = qx[i];
			ids[k] = qid[i];
			mat[k] = qmat[i];
			ver[k] = qver[i];
		}
		qa = a;
		qaBase = aBase;
		qm = mm;
		qb = b;
		qbBase = bBase;
		qn = nn;
		qrev = rev;
		qx = x;
		qid = ids;
		qmat = mat;
		qver = ver;
		head = 0;
	}

	/**
	 * Advances every lane by one row, reporting the extensions that end to {@code done} (which may submit more).
	 * Returns false, doing nothing, when no lane has an extension and none is queued.
	 */
	public boolean step(Done done) {
		t++;
		// Row setup: the next row of each lane, a new extension where one has passed its last row or where the lane
		// is idle.
		IntVector aiV = IntVector.fromArray(S, ai, 0).add(1);
		aiV.intoArray(ai, 0);
		if (aiV.compare(VectorOperators.GT, IntVector.fromArray(S, m, 0)).anyTrue()
				|| size > 0 && IntVector.fromArray(S, busy, 0).compare(VectorOperators.EQ, 0).anyTrue()) {
			for (int l = 0; l < L; l++) {
				if (busy[l] != 0 && ai[l] > m[l]) {
					finish(l, done);
				}
				if (busy[l] == 0 && size > 0) {
					start(l, t - 1);
					ai[l] = 1;
				}
			}
			aiV = IntVector.fromArray(S, ai, 0);
		}
		final IntVector bSizeV = IntVector.fromArray(S, bSize, 0);
		// The band's letters: cell bi reads B's letter bi - 1. Written ahead in chunks.
		if (IntVector.fromArray(S, written, 0).compare(VectorOperators.LT, bSizeV.sub(2)).anyTrue()) {
			for (int l = 0; l < L; l++) {
				if (written[l] < bSize[l] - 2) {
					letters(l, bSize[l] + 64);
				}
			}
		}
		final IntVector firstV = IntVector.fromArray(S, first, 0);
		final IntVector kStartV = firstV.sub(aiV).add(off), widthV = bSizeV.sub(firstV);
		final VectorMask<Integer> working = widthV.compare(VectorOperators.GT, 0);
		if (!working.anyTrue()) {
			return size > 0;
		}
		final int kLo = kStartV.reduceLanes(VectorOperators.MIN, working);
		final int kHi = kStartV.add(widthV).reduceLanes(VectorOperators.MAX, working);
		kStartV.intoArray(kStart, 0);
		widthV.intoArray(width, 0);
		// Room for the bands, the cells the gap may add after them, and the letters written ahead.
		if (kHi - kLo + reach + 72 > slots) {
			widen(kHi - kLo + reach + 72);
		}
		// The DP row's slots: the bands and the cells after them, and those of a new extension's first row.
		window(Math.min(kLo, off) - 2, Math.max(kHi, off) + reach + 3);
		IntVector.fromArray(S, aRing, (t & (aSize - 1)) * L).intoArray(row, 0);
		row(kLo, kHi);
		// Row end, as semiGapped: the next row's band ends after the last live cell, or after the cells that the
		// gap along the row still reaches, then gets a dead cell.
		final VectorMask<Integer> busyV = IntVector.fromArray(S, busy, 0).compare(VectorOperators.NE, 0);
		final VectorMask<Integer> ended = IntVector.fromArray(S, first, 0).compare(VectorOperators.EQ, bSizeV)
				.and(busyV);
		final VectorMask<Integer> going = busyV.andNot(ended);
		final IntVector nV = IntVector.fromArray(S, n, 0);
		final IntVector lastV = IntVector.fromArray(S, last, 0);
		final VectorMask<Integer> shrink = lastV.compare(VectorOperators.LT, bSizeV.sub(1));
		IntVector bs = bSizeV.blend(lastV.add(1), shrink);
		final IntVector floorV = IntVector.fromArray(S, best, 0).sub(IntVector.fromArray(S, xDrop, 0));
		if (going.andNot(shrink).and(IntVector.fromArray(S, gapRow, 0).compare(VectorOperators.GE, floorV))
				.and(bs.compare(VectorOperators.LE, nV)).anyTrue()) {
			bs.intoArray(tmp, 0);
			for (int l = 0; l < L; l++) {
				if (going.laneIsSet(l) && !shrink.laneIsSet(l)) {
					final int floor = best[l] - xDrop[l], a = ai[l], nl = n[l];
					int b = tmp[l];
					final byte[] tr = trace[t & (steps - 1)];
					final int lo = traceLo[t & (steps - 1)];
					for (int g = gapRow[l]; g >= floor && b <= nl; g -= ext, b++) {
						final int k = (b - a + off - dpBase) * L + l;
						dpBest[k] = g;
						dpGap[k] = g - goe;
						tr[(b - a + off - lo) * L + l] = (byte) DEL;
					}
					tmp[l] = b;
				}
			}
			bs = IntVector.fromArray(S, tmp, 0);
		}
		final VectorMask<Integer> dead = bs.compare(VectorOperators.LE, nV).and(going);
		bs.sub(aiV).add(off - dpBase).mul(L).add(IntVector.zero(S).addIndex(1)).intoArray(tmp, 0);
		for (int l = 0; l < L; l++) {
			if (dead.laneIsSet(l)) {
				dpBest[tmp[l]] = MININT;
				dpGap[tmp[l]] = MININT;
			}
		}
		bs.add(1, dead).blend(bSizeV, going.not()).intoArray(bSize, 0);
		if (ended.anyTrue()) {
			for (int l = 0; l < L; l++) {
				if (ended.laneIsSet(l)) {
					finish(l, done);
					if (size > 0) {
						start(l, t);
					}
				}
			}
		}
		return true;
	}

	/** The cells of slots [kLo, kHi) of the current rows of all lanes. */
	private void row(int kLo, int kHi) {
		final IntVector kStartV = IntVector.fromArray(S, kStart, 0), widthV = IntVector.fromArray(S, width, 0);
		final IntVector rowV = IntVector.fromArray(S, row, 0).add(LANE_BASE);
		final IntVector xDropV = IntVector.fromArray(S, xDrop, 0);
		final IntVector bestBefore = IntVector.fromArray(S, best, 0);
		final int goe = this.goe, ext = this.ext;
		IntVector bestV = bestBefore;
		IntVector bestK = IntVector.zero(S), firstK = IntVector.broadcast(S, Integer.MAX_VALUE),
				lastK = IntVector.broadcast(S, Integer.MIN_VALUE), gapRowV = IntVector.broadcast(S, MININT);
		final int[] db = dpBest, dg = dpGap, rg = ring, ix = idx;
		final int base = ((kLo + t) & (slots - 1)) - kLo, dp = dpBase;
		// This step's traceback bytes: the slots of the bands and those the row's gap may add after them.
		final int need = (kHi - kLo + reach + 4) * L;
		byte[] tr = trace[t & (steps - 1)];
		if (tr == null || tr.length < need) {
			// With room to spare, so that a step's array is rarely replaced.
			tr = new byte[Math.max(2 * need, 256 * L)];
			trace[t & (steps - 1)] = tr;
		}
		traceLo[t & (steps - 1)] = kLo;
		for (int k = kLo; k < kHi; k++) {
			final IntVector kv = IntVector.broadcast(S, k);
			// The cell's place in its lane's band: outside it the cell is dead, at 0 it has no diagonal (semiGapped
			// starts a row with MININT).
			final IntVector d = kv.sub(kStartV);
			rowV.or(IntVector.fromArray(S, rg, (base + k) * L)).intoArray(ix, 0);
			final IntVector diag = IntVector.fromArray(S, mats, 0, ix, 0);
			final int at = (k - dp) * L;
			final IntVector oldBest = IntVector.fromArray(S, db, at);
			final IntVector oldGap = IntVector.fromArray(S, dg, at + L);
			final IntVector withDiag = oldBest.add(diag).blend(MININT, d.compare(VectorOperators.EQ, 0));
			final IntVector s = withDiag.max(oldGap.max(gapRowV));
			// alignEx's choice (ties keep the earlier move): the gap along the column over the diagonal, then the
			// gap along the row over both.
			IntVector script = IntVector.broadcast(S, SUB).blend(INS, oldGap.compare(VectorOperators.GT, withDiag))
					.blend(DEL, gapRowV.compare(VectorOperators.GT, withDiag.max(oldGap)));
			// A cell above the best score is live (the floor is below the best), so the best is a max.
			final VectorMask<Integer> dead = s.compare(VectorOperators.LT, bestV.sub(xDropV))
					.or(d.compare(VectorOperators.UGE, widthV));
			bestK = bestK.blend(kv, s.compare(VectorOperators.GT, bestV).andNot(dead));
			bestV = bestV.max(s.blend(MININT, dead));
			firstK = firstK.min(kv.blend(Integer.MAX_VALUE, dead));
			lastK = lastK.blend(kv, dead.not());
			// A dead cell keeps the gap along the row, and its column's gap.
			final IntVector open = s.sub(goe);
			// A live cell's gaps that extend rather than open.
			script = script.add(EXTEND_GAP_B, oldGap.sub(ext).compare(VectorOperators.GE, open).andNot(dead))
					.add(EXTEND_GAP_A, gapRowV.sub(ext).compare(VectorOperators.GE, open).andNot(dead));
			script.convertShape(VectorOperators.I2B, B, 0).reinterpretAsBytes().intoArray(tr, (k - kLo) * L);
			gapRowV = open.max(gapRowV.sub(ext)).blend(gapRowV, dead);
			s.blend(MININT, dead).intoArray(db, at);
			open.max(oldGap.sub(ext)).blend(oldGap, dead).intoArray(dg, at);
		}
		// Back to semiGapped's state: bi = k + ai - off.
		final IntVector aiV = IntVector.fromArray(S, ai, 0);
		final IntVector toBi = aiV.sub(off);
		final VectorMask<Integer> improved = bestV.compare(VectorOperators.GT, bestBefore);
		bestV.intoArray(best, 0);
		IntVector.fromArray(S, aOffset, 0).blend(aiV, improved).intoArray(aOffset, 0);
		IntVector.fromArray(S, bOffset, 0).blend(bestK.add(toBi), improved).intoArray(bOffset, 0);
		// No live cell: first reaches bSize, which ends the extension.
		firstK.add(toBi)
				.blend(IntVector.fromArray(S, bSize, 0), firstK.compare(VectorOperators.EQ, Integer.MAX_VALUE))
				.intoArray(first, 0);
		lastK.add(toBi).intoArray(last, 0);
		gapRowV.intoArray(gapRow, 0);
	}

	/**
	 * Makes the DP row's window hold slots [lo, hi], moving it, and growing it where needed, with room to spare.
	 * Copies the slots in use of each lane, which lie between rows within its band (in the coordinates of its last or
	 * its next row) and the cells a row's gap adds after it.
	 */
	private void window(int lo, int hi) {
		if (lo >= dpBase && hi < dpBase + dpSlots) {
			return;
		}
		int n = dpSlots;
		while (n < 2 * (hi - lo + 1)) {
			n *= 2;
		}
		final int base = (lo + hi) / 2 - n / 2;
		final int[] best = new int[n * L], gap = new int[n * L];
		for (int l = 0; l < L; l++) {
			if (busy[l] != 0) {
				final int from = Math.max(first[l] - ai[l] + off - 3, Math.max(dpBase, base));
				final int to = Math.min(bSize[l] - ai[l] + off + reach + 3,
						Math.min(dpBase + dpSlots, base + n) - 1);
				for (int k = from; k <= to; k++) {
					best[(k - base) * L + l] = dpBest[(k - dpBase) * L + l];
					gap[(k - base) * L + l] = dpGap[(k - dpBase) * L + l];
				}
			}
		}
		dpBest = best;
		dpGap = gap;
		dpBase = base;
		dpSlots = n;
	}

	/**
	 * Makes the letter ring at least {@code needed} positions, and writes each lane's letters again from its band's
	 * start.
	 */
	private void widen(int needed) {
		int s = slots;
		while (s < needed) {
			s *= 2;
		}
		ring = new int[2 * s * L];
		slots = s;
		for (int l = 0; l < L; l++) {
			if (busy[l] != 0) {
				final int upTo = written[l];
				written[l] = Math.max(-1, first[l] - 3);
				letters(l, upTo);
			}
		}
	}

	/** Writes lane l's B letters up to letter {@code upTo} (as semiGapped reads them) into the ring. */
	private void letters(int l, int upTo) {
		final byte[] b = subject[l];
		final int nn = n[l], bb = bBase[l], position = shift[l] + off + 1;
		// As GappedAligner.bLetters: forwards b[bBase + 1 ..], backwards b[bBase + n - 1 ..], cut at b's ends.
		final int end = Math.min(upTo + 1,
				reverse[l] ? Math.min(nn + 1, bb + nn) : Math.min(nn + 1, b.length - bb - 1));
		for (int k = written[l] + 1; k < end; k++) {
			final int letter = (reverse[l] ? b[bb + nn - 1 - k] : b[bb + 1 + k]) & 31;
			final int p = (k + position) & (slots - 1);
			ring[p * L + l] = letter;
			ring[(p + slots) * L + l] = letter;
		}
		written[l] = Math.max(written[l], end - 1);
	}

	private void finish(int l, Done done) {
		final int i = id[l], score = best[l], a = aOffset[l], b = bOffset[l];
		walk(l, a, b);
		idle(l);
		done.done(i, score, a, b, ops);
	}

	/** alignEx's walk from the best cell (a, b) of lane l's extension back to the start, into {@link #ops}. */
	private void walk(int l, int a, int b) {
		ops.clear();
		final int sh = shift[l];
		int script = SUB;
		while (a > 0 || b > 0) {
			int next;
			if (a == 0) {
				next = DEL; // the first row
			} else {
				final int st = (a + sh) & (steps - 1);
				next = trace[st][(b - a + off - traceLo[st]) * L + l];
			}
			switch (script) {
			case DEL:
				script = next & 0x07;
				if ((next & EXTEND_GAP_A) != 0) {
					script = DEL;
				}
				break;
			case INS:
				script = next & 0x07;
				if ((next & EXTEND_GAP_B) != 0) {
					script = INS;
				}
				break;
			default:
				script = next & 0x07;
				break;
			}
			if (script == DEL) {
				b--;
			} else if (script == INS) {
				a--;
			} else {
				a--;
				b--;
			}
			ops.add(script, 1);
		}
	}

	/** A lane without an extension, which no step changes. */
	private void idle(int l) {
		busy[l] = 0;
		m[l] = Integer.MAX_VALUE;
		n[l] = -1;
		first[l] = 0;
		bSize[l] = 0;
		written[l] = Integer.MAX_VALUE;
		subject[l] = null;
	}

	/** Starts the next queued extension in lane l, whose row ai comes at step {@code ai + sh}. */
	private void start(int l, int sh) {
		final int i = head;
		head = (head + 1) % qa.length;
		size--;
		final byte[] a = qa[i];
		final int aBase = qaBase[i], mm = qm[i], nn = qn[i];
		final boolean rev = qrev[i];
		final int xd = Math.max(qx[i], goe);
		if (xd / ext + 2 + 72 > slots) {
			widen(xd / ext + 2 + 72);
		}
		window(off - 1, off + xd / ext + 3);
		qa[i] = null;
		subject[l] = qb[i];
		qb[i] = null;
		bBase[l] = qbBase[i];
		reverse[l] = rev;
		shift[l] = sh;
		written[l] = -1;
		for (int k = 1; k <= mm; k++) {
			aRing[((k + sh) & (aSize - 1)) * L + l] = ((rev ? a[aBase + mm - k] : a[aBase + k]) & 31) << 5;
		}
		// semiGapped's first row (ai = 0), in slots bi + off. Its band, like any, is less than the X-drop reach.
		int score = -goe;
		dpBest[(off - dpBase) * L + l] = 0;
		dpGap[(off - dpBase) * L + l] = -goe;
		int bi;
		for (bi = 1; bi <= nn; bi++) {
			if (score < -xd) {
				break;
			}
			final int k = (bi + off - dpBase) * L + l;
			dpBest[k] = score;
			dpGap[k] = score - goe;
			score -= ext;
		}
		if (laneMatrix[l] != qver[i]) {
			laneMatrix[l] = qver[i];
			System.arraycopy(qmat[i], 0, mats, l << 10, 1 << 10);
		}
		qmat[i] = null;
		id[l] = qid[i];
		ai[l] = 0;
		m[l] = mm;
		n[l] = nn;
		first[l] = 0;
		bSize[l] = bi;
		best[l] = 0;
		aOffset[l] = 0;
		bOffset[l] = 0;
		xDrop[l] = xd;
		busy[l] = -1;
	}
}
