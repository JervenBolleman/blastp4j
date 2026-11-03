package swiss.sib.blastp4j.align;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Many independent {@link GappedAligner#semiGapped} extensions (BLOSUM62) at once, one per lane of a vector of 16
 * shorts.
 * Each lane computes exactly what semiGapped computes for its extension: the same cells, the same X-drop, the same
 * band, the same best score and end point. Extensions are {@linkplain #submit submitted} to a queue; {@link #step}
 * advances every lane by one row, and a lane whose extension ends reports it and takes the next one from the queue.
 * <p>
 * The lanes share their columns by keeping the DP row in diagonal coordinates: lane l's cell (ai, bi) is in slot
 * {@code k = bi - ai + off}, at {@code k * L + l}. A band follows its diagonal, so the bands of all lanes mostly cover
 * the same slots, and a step over the slots of all bands handles one cell of each lane at once. In these coordinates
 * the cell above (ai - 1, bi) is in slot k + 1 and the diagonal one (ai - 1, bi - 1) in slot k, so a row is computed
 * in place, in increasing k. B's letters are kept in a ring indexed by slot plus row step, so that the letters of all
 * lanes for a slot are one load; only the score matrix is read with a gather. The DP row is a window on the slots,
 * which moves and grows with the bands, and the letter ring grows with their span, so that the buffers stay small:
 * arrays sized for the longest subject are humongous objects for G1, which then grows the heap several times over.
 * <p>
 * Scores and slots are shorts in the row loop: an extension whose A and B could make them overflow ({@link #fits})
 * is computed with semiGapped instead, and reported at the next step. The bookkeeping between rows is on ints, in two
 * vectors of 8.
 * <p>
 * Pays where the machine has 8 int lanes with a gather instruction (x86 with AVX2): {@link #AVAILABLE}. Not thread
 * safe: use one per thread.
 */
public final class SemiGappedLanes {
	/** Receives the result of an extension: semiGapped's score and its {@code aOffset} and {@code bOffset}. */
	public interface Done {
		void done(int id, int score, int aOffset, int bOffset);
	}

	/** The bookkeeping: half the lanes per vector. */
	private static final VectorSpecies<Integer> S = IntVector.SPECIES_256;
	/** The row loop: all lanes in one vector. */
	private static final VectorSpecies<Short> SS = ShortVector.SPECIES_256;
	/** Lanes. */
	public static final int L = SS.length();
	/** Lanes per bookkeeping vector. */
	private static final int H = S.length();
	/** The dead cell's score in the row loop: below any score, and far enough from the short's minimum. */
	private static final short MINSHORT = -16384;
	/** The longest A and the most letters of A and B together that the shorts allow ({@link #fits}). */
	private static final int SHORT_MAX_A = 1400, SHORT_MAX_AB = 30000;
	/**
	 * Whether to use the lanes: where the preferred vector has at least 8 ints (here a gather is an instruction), not
	 * in a native image (where Graal compiles a gather to scalar loads), and unless {@code -Dblastp4j.gappedLanes=false}.
	 */
	public static final boolean AVAILABLE = IntVector.SPECIES_PREFERRED.length() >= 8
			&& System.getProperty("org.graalvm.nativeimage.imagecode") == null
			&& Boolean.parseBoolean(System.getProperty("blastp4j.gappedLanes", "true"));
	private static final int[] MAT = GappedAligner.Rows.BLOSUM62;

	private final int goe, ext;
	/** The longest A an extension may have. */
	private final int maxA;
	/** The slot of cell (0, 0): diagonal coordinates are {@code bi - ai + off}. */
	private final int off;
	/** The DP row: slots {@code [dpBase, dpBase + dpSlots)}, slot k of lane l at {@code (k - dpBase) * L + l}. */
	private short[] dpBest, dpGap;
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
	private final short[] stmp = new short[L];
	/** Extensions done with semiGapped ({@link #fits}), reported at the next step: id, score, aOffset, bOffset. */
	private int[] doneNow = new int[16];
	private int doneCount;
	private final GappedAligner scalar;

	// The queue of submitted extensions.
	private byte[][] qa = new byte[16][], qb = new byte[16][];
	private int[] qaBase = new int[16], qm = new int[16], qbBase = new int[16], qn = new int[16], qx = new int[16],
			qid = new int[16];
	private boolean[] qrev = new boolean[16];
	private int head, size;

	/** Lanes for the given gap costs, and for extensions with A of at most {@code maxA} letters. */
	public SemiGappedLanes(int gapOpen, int gapExtend, int maxA) {
		this(gapOpen, gapExtend, maxA, 256);
	}

	/** With DP row and letter ring of {@code slots} slots to start with (the tests start small, to move and grow). */
	SemiGappedLanes(int gapOpen, int gapExtend, int maxA, int slots) {
		this.slots = slots;
		this.dpSlots = slots;
		goe = gapOpen + gapExtend;
		ext = gapExtend;
		this.maxA = maxA;
		off = maxA + 2;
		dpBest = new short[dpSlots * L];
		dpGap = new short[dpSlots * L];
		scalar = new GappedAligner(gapOpen, gapExtend);
		dpBase = off - dpSlots / 2;
		ring = new int[2 * slots * L];
		aSize = Integer.highestOneBit(maxA + 2) * 2;
		aRing = new int[aSize * L];
		for (int l = 0; l < L; l++) {
			idle(l);
		}
	}

	/** Extensions submitted and not yet in a lane or reported. */
	public int queued() {
		return size + doneCount / 4;
	}

	/**
	 * Whether an extension with A of m letters and B of n fits the row loop's shorts: its scores (at most 11 per letter
	 * of the shorter of A and B, BLOSUM62's best) and its slots (within m + n of each other) stay well within them.
	 */
	static boolean fits(int m, int n) {
		return m <= SHORT_MAX_A && m + n <= SHORT_MAX_AB;
	}

	/**
	 * Submits {@code semiGapped(a, aBase, m, b, bBase, n, reverse, xDrop)}, with {@code 0 < m <= maxA} and
	 * {@code n > 0}; {@link #step} reports its result under {@code id}.
	 */
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id) {
		if (m <= 0 || n <= 0 || m > maxA) {
			throw new IllegalArgumentException("m " + m + ", n " + n);
		}
		if (!fits(m, n)) {
			final int score = scalar.semiGapped(a, aBase, m, b, bBase, n, reverse, xDrop);
			if (doneCount == doneNow.length) {
				doneNow = java.util.Arrays.copyOf(doneNow, 2 * doneCount);
			}
			doneNow[doneCount++] = id;
			doneNow[doneCount++] = score;
			doneNow[doneCount++] = scalar.aOffset;
			doneNow[doneCount++] = scalar.bOffset;
			return;
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
		size++;
	}

	private void grow() {
		final int c = qa.length, c2 = 2 * c;
		final byte[][] a = new byte[c2][], b = new byte[c2][];
		final int[] aBase = new int[c2], mm = new int[c2], bBase = new int[c2], nn = new int[c2], x = new int[c2],
				ids = new int[c2];
		final boolean[] rev = new boolean[c2];
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
		head = 0;
	}

	/**
	 * Advances every lane by one row, reporting the extensions that end to {@code done} (which may submit more).
	 * Returns false, doing nothing, when no lane has an extension and none is queued.
	 */
	public boolean step(Done done) {
		if (doneCount > 0) {
			final int[] d = doneNow;
			final int count = doneCount;
			doneNow = new int[Math.max(16, count)];
			doneCount = 0;
			for (int k = 0; k < count; k += 4) {
				done.done(d[k], d[k + 1], d[k + 2], d[k + 3]);
			}
			return true;
		}
		t++;
		// Row setup: the next row of each lane, a new extension where one has passed its last row or where the lane
		// is idle.
		boolean setup = false;
		for (int h = 0; h < L; h += H) {
			final IntVector aiV = IntVector.fromArray(S, ai, h).add(1);
			aiV.intoArray(ai, h);
			setup |= aiV.compare(VectorOperators.GT, IntVector.fromArray(S, m, h)).anyTrue()
					|| size > 0 && IntVector.fromArray(S, busy, h).compare(VectorOperators.EQ, 0).anyTrue();
		}
		if (setup) {
			for (int l = 0; l < L; l++) {
				if (busy[l] != 0 && ai[l] > m[l]) {
					finish(l, done);
				}
				if (busy[l] == 0 && size > 0) {
					start(l, t - 1);
					ai[l] = 1;
				}
			}
		}
		// The band's letters: cell bi reads B's letter bi - 1. Written ahead in chunks.
		boolean write = false;
		for (int h = 0; h < L; h += H) {
			write |= IntVector.fromArray(S, written, h)
					.compare(VectorOperators.LT, IntVector.fromArray(S, bSize, h).sub(2)).anyTrue();
		}
		if (write) {
			for (int l = 0; l < L; l++) {
				if (written[l] < bSize[l] - 2) {
					letters(l, bSize[l] + 64);
				}
			}
		}
		int kLo = Integer.MAX_VALUE, kHi = Integer.MIN_VALUE;
		for (int h = 0; h < L; h += H) {
			final IntVector firstV = IntVector.fromArray(S, first, h);
			final IntVector kStartV = firstV.sub(IntVector.fromArray(S, ai, h)).add(off);
			final IntVector widthV = IntVector.fromArray(S, bSize, h).sub(firstV);
			final VectorMask<Integer> working = widthV.compare(VectorOperators.GT, 0);
			kLo = Math.min(kLo, kStartV.reduceLanes(VectorOperators.MIN, working));
			kHi = Math.max(kHi, kStartV.add(widthV).reduceLanes(VectorOperators.MAX, working));
			kStartV.intoArray(kStart, h);
			widthV.intoArray(width, h);
		}
		if (kLo > kHi) {
			return size > 0 || doneCount > 0;
		}
		// Room for the bands, the cells the gap may add after them, and the letters written ahead.
		if (kHi - kLo + reach + 72 > slots) {
			widen(kHi - kLo + reach + 72);
		}
		// The DP row's slots: the bands and the cells after them, and those of a new extension's first row.
		window(Math.min(kLo, off) - 2, Math.max(kHi, off) + reach + 3);
		System.arraycopy(aRing, (t & (aSize - 1)) * L, row, 0, L);
		row(kLo, kHi);
		// Row end, as semiGapped: the next row's band ends after the last live cell, or after the cells that the
		// gap along the row still reaches, then gets a dead cell.
		boolean anyEnded = false;
		for (int h = 0; h < L; h += H) {
			final IntVector bSizeV = IntVector.fromArray(S, bSize, h);
			final IntVector aiV = IntVector.fromArray(S, ai, h);
			final VectorMask<Integer> busyV = IntVector.fromArray(S, busy, h).compare(VectorOperators.NE, 0);
			final VectorMask<Integer> ended = IntVector.fromArray(S, first, h).compare(VectorOperators.EQ, bSizeV)
					.and(busyV);
			anyEnded |= ended.anyTrue();
			final VectorMask<Integer> going = busyV.andNot(ended);
			final IntVector nV = IntVector.fromArray(S, n, h);
			final IntVector lastV = IntVector.fromArray(S, last, h);
			final VectorMask<Integer> shrink = lastV.compare(VectorOperators.LT, bSizeV.sub(1));
			IntVector bs = bSizeV.blend(lastV.add(1), shrink);
			final IntVector floorV = IntVector.fromArray(S, best, h).sub(IntVector.fromArray(S, xDrop, h));
			if (going.andNot(shrink).and(IntVector.fromArray(S, gapRow, h).compare(VectorOperators.GE, floorV))
					.and(bs.compare(VectorOperators.LE, nV)).anyTrue()) {
				bs.intoArray(tmp, h);
				for (int l = h; l < h + H; l++) {
					if (going.laneIsSet(l - h) && !shrink.laneIsSet(l - h)) {
						final int floor = best[l] - xDrop[l], a = ai[l], nl = n[l];
						int b = tmp[l];
						for (int g = gapRow[l]; g >= floor && b <= nl; g -= ext, b++) {
							final int k = (b - a + off - dpBase) * L + l;
							dpBest[k] = (short) g;
							dpGap[k] = (short) (g - goe);
						}
						tmp[l] = b;
					}
				}
				bs = IntVector.fromArray(S, tmp, h);
			}
			final VectorMask<Integer> dead = bs.compare(VectorOperators.LE, nV).and(going);
			bs.sub(aiV).add(off - dpBase).mul(L).add(IntVector.zero(S).addIndex(1).add(h)).intoArray(tmp, h);
			for (int l = h; l < h + H; l++) {
				if (dead.laneIsSet(l - h)) {
					dpBest[tmp[l]] = MINSHORT;
					dpGap[tmp[l]] = MINSHORT;
				}
			}
			bs.add(1, dead).blend(bSizeV, going.not()).intoArray(bSize, h);
		}
		if (anyEnded) {
			for (int l = 0; l < L; l++) {
				if (busy[l] != 0 && first[l] == bSize[l]) {
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
		// In the loop, slots count from kLo.
		final ShortVector kStartV = shorts(kStart, -kLo), widthV = shorts(width, 0);
		final ShortVector xDropV = shorts(xDrop, 0);
		final ShortVector bestBefore = shorts(best, 0);
		final short goe = (short) this.goe, ext = (short) this.ext;
		ShortVector bestV = bestBefore;
		ShortVector bestK = ShortVector.zero(SS), firstK = ShortVector.broadcast(SS, Short.MAX_VALUE),
				lastK = ShortVector.broadcast(SS, Short.MIN_VALUE), gapRowV = ShortVector.broadcast(SS, MINSHORT);
		final short[] db = dpBest, dg = dpGap;
		final int[] rg = ring, ix = idx, rw = row;
		final int base = ((kLo + t) & (slots - 1)) - kLo, dp = dpBase;
		for (int k = kLo; k < kHi; k++) {
			final ShortVector kv = ShortVector.broadcast(SS, (short) (k - kLo));
			// The cell's place in its lane's band: outside it the cell is dead, at 0 it has no diagonal (semiGapped
			// starts a row with MININT).
			final ShortVector d = kv.sub(kStartV);
			final int r = (base + k) * L;
			IntVector.fromArray(S, rw, 0).or(IntVector.fromArray(S, rg, r)).intoArray(ix, 0);
			IntVector.fromArray(S, rw, H).or(IntVector.fromArray(S, rg, r + H)).intoArray(ix, H);
			final ShortVector score = (ShortVector) IntVector.fromArray(S, MAT, 0, ix, 0)
					.convertShape(VectorOperators.I2S, SS, 0)
					.lanewise(VectorOperators.OR,
							IntVector.fromArray(S, MAT, 0, ix, H).convertShape(VectorOperators.I2S, SS, -1));
			final ShortVector diag = score.blend(MINSHORT, d.compare(VectorOperators.EQ, 0));
			final int at = (k - dp) * L;
			final ShortVector oldBest = ShortVector.fromArray(SS, db, at);
			final ShortVector oldGap = ShortVector.fromArray(SS, dg, at + L);
			final ShortVector s = oldBest.add(diag).max(oldGap.max(gapRowV));
			// A cell above the best score is live (the floor is below the best), so the best is a max.
			final VectorMask<Short> dead = s.compare(VectorOperators.LT, bestV.sub(xDropV))
					.or(d.compare(VectorOperators.UGE, widthV));
			bestK = bestK.blend(kv, s.compare(VectorOperators.GT, bestV).andNot(dead));
			bestV = bestV.max(s.blend(MINSHORT, dead));
			firstK = firstK.min(kv.blend(Short.MAX_VALUE, dead));
			lastK = lastK.blend(kv, dead.not());
			// A dead cell keeps the gap along the row, and its column's gap.
			final ShortVector open = s.sub(goe);
			gapRowV = open.max(gapRowV.sub(ext)).blend(gapRowV, dead);
			s.blend(MINSHORT, dead).intoArray(db, at);
			open.max(oldGap.sub(ext)).blend(oldGap, dead).intoArray(dg, at);
		}
		// Back to semiGapped's state, in ints: bi = k + ai - off.
		final VectorMask<Short> improvedS = bestV.compare(VectorOperators.GT, bestBefore);
		final VectorMask<Short> noneLive = firstK.compare(VectorOperators.EQ, Short.MAX_VALUE);
		bestV.intoArray(stmp, 0);
		for (int l = 0; l < L; l++) {
			best[l] = stmp[l];
		}
		bestK.intoArray(stmp, 0);
		for (int l = 0; l < L; l++) {
			if (improvedS.laneIsSet(l)) {
				aOffset[l] = ai[l];
				bOffset[l] = stmp[l] + kLo + ai[l] - off;
			}
		}
		// No live cell: first reaches bSize, which ends the extension.
		firstK.intoArray(stmp, 0);
		for (int l = 0; l < L; l++) {
			first[l] = noneLive.laneIsSet(l) ? bSize[l] : stmp[l] + kLo + ai[l] - off;
		}
		lastK.intoArray(stmp, 0);
		for (int l = 0; l < L; l++) {
			last[l] = stmp[l] == Short.MIN_VALUE ? Integer.MIN_VALUE : stmp[l] + kLo + ai[l] - off;
		}
		gapRowV.intoArray(stmp, 0);
		for (int l = 0; l < L; l++) {
			gapRow[l] = stmp[l];
		}
	}

	/** The lanes' values plus {@code add}, as shorts. */
	private ShortVector shorts(int[] values, int add) {
		for (int l = 0; l < L; l++) {
			stmp[l] = (short) (values[l] + add);
		}
		return ShortVector.fromArray(SS, stmp, 0);
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
		final short[] best = new short[n * L], gap = new short[n * L];
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
		idle(l);
		done.done(i, score, a, b);
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
		dpGap[(off - dpBase) * L + l] = (short) -goe;
		int bi;
		for (bi = 1; bi <= nn; bi++) {
			if (score < -xd) {
				break;
			}
			final int k = (bi + off - dpBase) * L + l;
			dpBest[k] = (short) score;
			dpGap[k] = (short) (score - goe);
			score -= ext;
		}
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
