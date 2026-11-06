package swiss.sib.blastp4j.seed;

import java.util.List;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * The ungapped extensions of a batch's word hits, {@link #LANES} at a time ({@code s_BlastAaExtendTwoHit}, as
 * {@link WordFinder} does one at a time). An extension decides its diagonal's state, which later hits on that
 * diagonal read, so a hit is only queued here: its diagonal is marked {@link #PENDING}, and a later hit on it, a full
 * queue or the end of the subject {@link #flush}es the queue first. The queued hits are then on different diagonals,
 * independent of each other: their extensions run as lanes, one step of each lane per round, and their results are
 * applied in the order the hits came, as without the queue.
 */
final class UngappedLanes {
	/**
	 * A query's word finder, as the queue sees it: an interface, as {@link WordFinder}'s hidden copies are not
	 * {@link WordFinder}s.
	 */
	interface Owner {
		/** The query's rows of {@link WordFinder#PADDED}, by query offset. */
		int[] rows();

		int xDrop();

		/** The query's offset in the batch's concatenation. */
		int base();

		int queryLength();

		/** Applies the result of hit {@code (q, s)}'s extension, as the finder applies its own. */
		void extended(int q, int s, int diag, int score, int hq, int hs, int hl, boolean right, int lastSubject,
				List<UngappedHsp> out);
	}

	/** The machine's preferred vector of ints, at most 256 bits. */
	private static final VectorSpecies<Integer> I = IntVector.SPECIES_PREFERRED.length() > 8 ? IntVector.SPECIES_256
			: IntVector.SPECIES_PREFERRED;
	/** The extensions per flush: a lane each. */
	static final int LANES = I.length();
	/** {@link WordFinder#PADDED} as ints, for the gathers. */
	private static final int[] PADDED_INTS = new int[WordFinder.PADDED.length];
	static {
		for (int k = 0; k < PADDED_INTS.length; k++) {
			PADDED_INTS[k] = WordFinder.PADDED[k];
		}
	}
	/** The diagonal state of a queued hit: never a position (those stay below bit 29), and {@code EXTENDED} set. */
	static final int PENDING = Integer.MAX_VALUE;

	private final Owner[] finder = new Owner[LANES];
	@SuppressWarnings("unchecked")
	private final List<UngappedHsp>[] out = new List[LANES];
	private final int[] q = new int[LANES], s = new int[LANES], sLeft = new int[LANES], diag = new int[LANES];
	private int size;
	/** The rows of the batch's queries concatenated, as the finders' diagonals: each query from its base on. */
	private final int[] allRows;
	/** The current subject's letters as ints, for the gathers, and which subject they are. */
	private int[] subjectInts = new int[0];
	private byte[] subjectArray;
	private int subjectStart = -1;
	private final int[] qi = new int[LANES], si = new int[LANES];
	private final int[] qmap = new int[LANES], smap = new int[LANES], tmp = new int[LANES], xd = new int[LANES];

	UngappedLanes(int[] allRows) {
		this.allRows = allRows;
	}

	/** A new subject: the letters are converted again at the next flush. */
	void newSubject() {
		subjectArray = null;
	}

	// Per lane: the extension's start after the word, and the left and right extension's state.
	private final int[] qR = new int[LANES], sR = new int[LANES];
	private final int[] n = new int[LANES], score = new int[LANES], maxScore = new int[LANES], best = new int[LANES],
			i = new int[LANES];
	private final int[] leftScore = new int[LANES], leftD = new int[LANES], rightScore = new int[LANES],
			rightD = new int[LANES], last = new int[LANES];
	private final boolean[] right = new boolean[LANES];

	/** Queues hit {@code (q, s)} of {@code f} on diagonal {@code diag}; the first word ends at {@code sLeft}. */
	boolean add(Owner f, List<UngappedHsp> o, int qo, int so, int sl, int d) {
		finder[size] = f;
		out[size] = o;
		q[size] = qo;
		s[size] = so;
		sLeft[size] = sl;
		diag[size] = d;
		return ++size == LANES;
	}

	/** Extends the queued hits of the subject {@code residues[start .. start + length)} and applies the results. */
	void flush(byte[] residues, int start, int length) {
		final int m = size;
		if (m == 0) {
			return;
		}
		if (subjectArray != residues || subjectStart != start) {
			if (subjectInts.length < length) {
				subjectInts = new int[Math.max(length, 2 * subjectInts.length)];
			}
			for (int k = 0; k < length; k++) {
				subjectInts[k] = residues[start + k];
			}
			subjectArray = residues;
			subjectStart = start;
		}
		// The word's best starting score, then the left extensions: s_BlastAaExtendLeft with a start score of 0.
		for (int l = 0; l < m; l++) {
			final int[] rows = finder[l].rows();
			int sc = 0, ls = 0, d = 0;
			for (int k = 0; k < LookupTable.WORD_SIZE; k++) {
				sc += WordFinder.PADDED[(rows[q[l] + k] + residues[start + s[l] + k]) & 1023];
				if (sc > ls) {
					ls = sc;
					d = k + 1;
				}
			}
			qR[l] = q[l] + d;
			sR[l] = s[l] + d;
			n[l] = Math.min(sR[l] - 1, qR[l] - 1);
			xd[l] = finder[l].xDrop();
			qmap[l] = finder[l].base() + qR[l] - 1;
			smap[l] = sR[l] - 1;
		}
		for (int l = m; l < LANES; l++) {
			qmap[l] = 0;
			smap[l] = 0;
		}
		final IntVector xDrops = IntVector.fromArray(I, xd, 0);
		VectorMask<Integer> act = IntVector.fromArray(I, n, 0).compare(VectorOperators.GE, 0)
				.and(I.indexInRange(0, m));
		IntVector sc = IntVector.zero(I), max = IntVector.zero(I), fl = xDrops.neg();
		IntVector bestK = IntVector.broadcast(I, -1), nv = IntVector.fromArray(I, n, 0);
		final IntVector qv = IntVector.fromArray(I, qmap, 0), sv = IntVector.fromArray(I, smap, 0);
		for (int k = 0; act.anyTrue(); k++) {
			// Lane l reads query position qmap[l] - k and subject position smap[l] - k.
			// A gather checks the indices of all lanes, masked or not: those of finished lanes are kept in bounds.
			qv.sub(k).max(0).intoArray(qi, 0);
			sv.sub(k).max(0).intoArray(si, 0);
			final IntVector rows = IntVector.fromArray(I, allRows, 0, qi, 0);
			final IntVector res = IntVector.fromArray(I, subjectInts, 0, si, 0);
			rows.add(res).and(1023).intoArray(tmp, 0);
			sc = sc.add(IntVector.fromArray(I, PADDED_INTS, 0, tmp, 0, act));
			final VectorMask<Integer> gt = sc.compare(VectorOperators.GT, max).and(act);
			max = max.blend(sc, gt);
			bestK = bestK.blend(k, gt);
			fl = fl.blend(sc.sub(xDrops), gt);
			act = act.andNot(sc.compare(VectorOperators.LE, fl)).and(nv.compare(VectorOperators.GT, k));
		}
		max.intoArray(maxScore, 0);
		// best = n - k of the best round, n + 1 when no round improved.
		nv.sub(bestK).blend(nv.add(1), bestK.compare(VectorOperators.LT, 0)).intoArray(best, 0);
		// The right extensions, s_BlastAaExtendRight, of the hits whose left extension reached the first word.
		for (int l = 0; l < m; l++) {
			leftScore[l] = maxScore[l];
			leftD[l] = n[l] - best[l] + 1;
			right[l] = leftD[l] >= sR[l] - sLeft[l];
			rightScore[l] = 0;
			rightD[l] = 0;
			last[l] = sR[l];
			if (right[l]) {
				n[l] = Math.min(length - sR[l], finder[l].queryLength() - qR[l]);
				score[l] = maxScore[l];
			} else {
				n[l] = 0;
			}
		}
		for (int l = 0; l < m; l++) {
			qmap[l] = Math.min(finder[l].base() + qR[l], allRows.length - 1);
			smap[l] = Math.min(sR[l], subjectInts.length - 1);
		}
		VectorMask<Integer> ract = IntVector.fromArray(I, n, 0).compare(VectorOperators.GT, 0)
				.and(I.indexInRange(0, m));
		IntVector rsc = IntVector.fromArray(I, score, 0), rmax = rsc;
		IntVector rfl = rmax.sub(xDrops).max(0), rbest = IntVector.broadcast(I, -1);
		final IntVector rn = IntVector.fromArray(I, n, 0);
		IntVector at = rn;
		final IntVector rqv = IntVector.fromArray(I, qmap, 0).max(0), rsv = IntVector.fromArray(I, smap, 0).max(0);
		for (int k = 0; ract.anyTrue(); k++) {
			rqv.add(k).min(allRows.length - 1).intoArray(qi, 0);
			rsv.add(k).min(subjectInts.length - 1).intoArray(si, 0);
			final IntVector rows = IntVector.fromArray(I, allRows, 0, qi, 0);
			final IntVector res = IntVector.fromArray(I, subjectInts, 0, si, 0);
			rows.add(res).and(1023).intoArray(tmp, 0);
			rsc = rsc.add(IntVector.fromArray(I, PADDED_INTS, 0, tmp, 0, ract));
			final VectorMask<Integer> gt = rsc.compare(VectorOperators.GT, rmax).and(ract);
			rmax = rmax.blend(rsc, gt);
			rbest = rbest.blend(k, gt);
			rfl = rfl.blend(rsc.sub(xDrops).max(0), gt);
			// The step of a break is the extension's end; a lane that runs to n ends at n.
			final VectorMask<Integer> brk = rsc.compare(VectorOperators.LE, rfl).and(ract);
			at = at.blend(k, brk);
			ract = ract.andNot(brk).and(rn.compare(VectorOperators.GT, k + 1));
		}
		rmax.intoArray(maxScore, 0);
		rbest.intoArray(best, 0);
		at.intoArray(i, 0);
		for (int l = 0; l < m; l++) {
			if (right[l]) {
				rightScore[l] = maxScore[l];
				rightD[l] = best[l] + 1;
				last[l] = sR[l] + i[l];
			}
			finder[l].extended(q[l], s[l], diag[l], Math.max(leftScore[l], rightScore[l]), qR[l] - leftD[l],
					sR[l] - leftD[l], leftD[l] + rightD[l], right[l], last[l], out[l]);
			finder[l] = null;
			out[l] = null;
		}
		size = 0;
	}
}
