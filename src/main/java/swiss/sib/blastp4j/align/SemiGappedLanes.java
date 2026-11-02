package swiss.sib.blastp4j.align;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Many independent {@link GappedAligner#semiGapped} extensions (BLOSUM62) at once, one per lane of a vector of 8 ints.
 * Each lane computes exactly what semiGapped computes for its extension: the same cells, the same X-drop, the same
 * band, the same best score and end point. Extensions are {@linkplain #submit submitted} to a queue; {@link #step}
 * advances every lane by one row, and a lane whose extension ends reports it and takes the next one from the queue.
 * <p>
 * The lanes share their columns by keeping the DP row in diagonal coordinates: lane l's cell (ai, bi) is in slot
 * {@code k = bi - ai + off}, at {@code k * L + l}. A band follows its diagonal, so the bands of all lanes mostly cover
 * the same slots, and a step over the slots of all bands handles one cell of each lane at once. In these coordinates
 * the cell above (ai - 1, bi) is in slot k + 1 and the diagonal one (ai - 1, bi - 1) in slot k, so a row is computed
 * in place, in increasing k. B's letters are kept in a ring indexed by slot plus row step, so that the letters of all
 * lanes for a slot are one load; only the score matrix is read with a gather.
 * <p>
 * Pays where the machine has 8 int lanes with a gather instruction (x86 with AVX2): {@link #AVAILABLE}. Not thread
 * safe: use one per thread.
 */
public final class SemiGappedLanes {
	/** Receives the result of an extension: semiGapped's score and its {@code aOffset} and {@code bOffset}. */
	public interface Done {
		void done(int id, int score, int aOffset, int bOffset);
	}

	private static final VectorSpecies<Integer> S = IntVector.SPECIES_256;
	/** Lanes. */
	public static final int L = S.length();
	/**
	 * Whether to use the lanes: where the preferred vector has at least 8 ints (here a gather is an instruction), not
	 * in a native image (where Graal compiles a gather to scalar loads), and unless {@code -Dblastp4j.gappedLanes=false}.
	 */
	public static final boolean AVAILABLE = IntVector.SPECIES_PREFERRED.length() >= 8
			&& System.getProperty("org.graalvm.nativeimage.imagecode") == null
			&& Boolean.parseBoolean(System.getProperty("blastp4j.gappedLanes", "true"));
	private static final int MININT = GappedAligner.MININT;
	private static final int[] MAT = GappedAligner.Rows.BLOSUM62;

	private final int goe, ext;
	/** The longest A and B an extension may have. */
	private final int maxA, maxB;
	/** The slot of cell (0, 0): diagonal coordinates are {@code bi - ai + off}. */
	private final int off;
	private final int[] dpBest, dpGap;
	/**
	 * B's letters by lane, in a ring of {@code ringSize} positions, twice, so that a run of up to {@code ringSize}
	 * positions from any start is contiguous: lane l's letter for slot k in row step t is at position {@code k + t}.
	 * An extension whose row ai comes at step {@code ai + shift} has letter bi (as semiGapped reads it) at position
	 * {@code bi + shift + off + 1}.
	 */
	private final int[] ring;
	/** A's letters as matrix rows ({@code letter << 5}) by lane: row ai at position {@code ai + shift}. */
	private final int[] aRing;
	private final int ringSize;
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

	// The queue of submitted extensions.
	private byte[][] qa = new byte[16][], qb = new byte[16][];
	private int[] qaBase = new int[16], qm = new int[16], qbBase = new int[16], qn = new int[16], qx = new int[16],
			qid = new int[16];
	private boolean[] qrev = new boolean[16];
	private int head, size;

	/**
	 * Lanes for the given gap costs, and for extensions with A of at most {@code maxA} letters and B of at most
	 * {@code maxB}.
	 */
	public SemiGappedLanes(int gapOpen, int gapExtend, int maxA, int maxB) {
		goe = gapOpen + gapExtend;
		ext = gapExtend;
		this.maxA = maxA;
		this.maxB = maxB;
		off = maxA + 2;
		final int slots = off + maxB + 4;
		dpBest = new int[slots * L];
		dpGap = new int[slots * L];
		ringSize = Integer.highestOneBit(maxA + maxB + 4) * 2;
		ring = new int[2 * ringSize * L];
		aRing = new int[ringSize * L];
		for (int l = 0; l < L; l++) {
			idle(l);
		}
	}

	/** Extensions submitted and not yet in a lane. */
	public int queued() {
		return size;
	}

	/**
	 * Submits {@code semiGapped(a, aBase, m, b, bBase, n, reverse, xDrop)}, with {@code 0 < m <= maxA} and
	 * {@code 0 < n <= maxB}; {@link #step} reports its result under {@code id}.
	 */
	public void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id) {
		if (m <= 0 || n <= 0 || m > maxA || n > maxB) {
			throw new IllegalArgumentException("m " + m + ", n " + n);
		}
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
		IntVector.fromArray(S, aRing, (t & (ringSize - 1)) * L).intoArray(row, 0);
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
					for (int g = gapRow[l]; g >= floor && b <= nl; g -= ext, b++) {
						final int k = b - a + off;
						dpBest[k * L + l] = g;
						dpGap[k * L + l] = g - goe;
					}
					tmp[l] = b;
				}
			}
			bs = IntVector.fromArray(S, tmp, 0);
		}
		final VectorMask<Integer> dead = bs.compare(VectorOperators.LE, nV).and(going);
		bs.sub(aiV).add(off).mul(L).add(IntVector.zero(S).addIndex(1)).intoArray(tmp, 0);
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
		final IntVector rowV = IntVector.fromArray(S, row, 0);
		final IntVector xDropV = IntVector.fromArray(S, xDrop, 0);
		final IntVector bestBefore = IntVector.fromArray(S, best, 0);
		final int goe = this.goe, ext = this.ext;
		IntVector bestV = bestBefore;
		IntVector bestK = IntVector.zero(S), firstK = IntVector.broadcast(S, Integer.MAX_VALUE),
				lastK = IntVector.broadcast(S, Integer.MIN_VALUE), gapRowV = IntVector.broadcast(S, MININT);
		final int[] db = dpBest, dg = dpGap, rg = ring, ix = idx;
		final int base = ((kLo + t) & (ringSize - 1)) - kLo;
		for (int k = kLo; k < kHi; k++) {
			final IntVector kv = IntVector.broadcast(S, k);
			// The cell's place in its lane's band: outside it the cell is dead, at 0 it has no diagonal (semiGapped
			// starts a row with MININT).
			final IntVector d = kv.sub(kStartV);
			rowV.or(IntVector.fromArray(S, rg, (base + k) * L)).intoArray(ix, 0);
			final IntVector diag = IntVector.fromArray(S, MAT, 0, ix, 0)
					.blend(MININT, d.compare(VectorOperators.EQ, 0));
			final IntVector oldBest = IntVector.fromArray(S, db, k * L);
			final IntVector oldGap = IntVector.fromArray(S, dg, (k + 1) * L);
			final IntVector s = oldBest.add(diag).max(oldGap.max(gapRowV));
			// A cell above the best score is live (the floor is below the best), so the best is a max.
			final VectorMask<Integer> dead = s.compare(VectorOperators.LT, bestV.sub(xDropV))
					.or(d.compare(VectorOperators.UGE, widthV));
			bestK = bestK.blend(kv, s.compare(VectorOperators.GT, bestV).andNot(dead));
			bestV = bestV.max(s.blend(MININT, dead));
			firstK = firstK.min(kv.blend(Integer.MAX_VALUE, dead));
			lastK = lastK.blend(kv, dead.not());
			// A dead cell keeps the gap along the row, and its column's gap.
			final IntVector open = s.sub(goe);
			gapRowV = open.max(gapRowV.sub(ext)).blend(gapRowV, dead);
			s.blend(MININT, dead).intoArray(db, k * L);
			open.max(oldGap.sub(ext)).blend(oldGap, dead).intoArray(dg, k * L);
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

	/** Writes lane l's B letters up to letter {@code upTo} (as semiGapped reads them) into the ring. */
	private void letters(int l, int upTo) {
		final byte[] b = subject[l];
		final int nn = n[l], bb = bBase[l], position = shift[l] + off + 1;
		// As GappedAligner.bLetters: forwards b[bBase + 1 ..], backwards b[bBase + n - 1 ..], cut at b's ends.
		final int end = Math.min(upTo + 1,
				reverse[l] ? Math.min(nn + 1, bb + nn) : Math.min(nn + 1, b.length - bb - 1));
		for (int k = written[l] + 1; k < end; k++) {
			final int letter = (reverse[l] ? b[bb + nn - 1 - k] : b[bb + 1 + k]) & 31;
			final int p = (k + position) & (ringSize - 1);
			ring[p * L + l] = letter;
			ring[(p + ringSize) * L + l] = letter;
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
		qa[i] = null;
		subject[l] = qb[i];
		qb[i] = null;
		bBase[l] = qbBase[i];
		reverse[l] = rev;
		shift[l] = sh;
		written[l] = -1;
		for (int k = 1; k <= mm; k++) {
			aRing[((k + sh) & (ringSize - 1)) * L + l] = ((rev ? a[aBase + mm - k] : a[aBase + k]) & 31) << 5;
		}
		// semiGapped's first row (ai = 0), in slots bi + off.
		int score = -goe;
		dpBest[off * L + l] = 0;
		dpGap[off * L + l] = -goe;
		int bi;
		for (bi = 1; bi <= nn; bi++) {
			if (score < -xd) {
				break;
			}
			dpBest[(bi + off) * L + l] = score;
			dpGap[(bi + off) * L + l] = score - goe;
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
