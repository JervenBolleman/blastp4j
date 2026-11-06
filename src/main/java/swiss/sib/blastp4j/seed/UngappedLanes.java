package swiss.sib.blastp4j.seed;

import java.util.List;

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

		int queryLength();

		/** Applies the result of hit {@code (q, s)}'s extension, as the finder applies its own. */
		void extended(int q, int s, int diag, int score, int hq, int hs, int hl, boolean right, int lastSubject,
				List<UngappedHsp> out);
	}

	/** The extensions per flush. */
	static final int LANES = 8;
	/** The diagonal state of a queued hit: never a position (those stay below bit 29), and {@code EXTENDED} set. */
	static final int PENDING = Integer.MAX_VALUE;

	private final Owner[] finder = new Owner[LANES];
	@SuppressWarnings("unchecked")
	private final List<UngappedHsp>[] out = new List[LANES];
	private final int[] q = new int[LANES], s = new int[LANES], sLeft = new int[LANES], diag = new int[LANES];
	private int size;

	// Per lane: the extension's start after the word, and the left and right extension's state.
	private final int[] qR = new int[LANES], sR = new int[LANES];
	private final int[] n = new int[LANES], score = new int[LANES], maxScore = new int[LANES], best = new int[LANES],
			floor = new int[LANES], i = new int[LANES];
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
		// The word's best starting score, then the left extensions: s_BlastAaExtendLeft with a start score of 0.
		int rounds = 0;
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
			best[l] = n[l] + 1;
			score[l] = 0;
			maxScore[l] = 0;
			floor[l] = -finder[l].xDrop();
			i[l] = n[l];
			rounds = Math.max(rounds, n[l] + 1);
		}
		for (int k = 0; k < rounds; k++) {
			boolean any = false;
			for (int l = 0; l < m; l++) {
				final int at = i[l];
				if (at >= 0) {
					final int sc = score[l] + WordFinder.PADDED[(finder[l].rows()[qR[l] - 1 - n[l] + at]
							+ residues[start + sR[l] - 1 - n[l] + at]) & 1023];
					score[l] = sc;
					if (sc > maxScore[l]) {
						maxScore[l] = sc;
						best[l] = at;
						floor[l] = sc - finder[l].xDrop();
					}
					i[l] = sc <= floor[l] ? -1 : at - 1;
					any = true;
				}
			}
			if (!any) {
				break;
			}
		}
		// The right extensions, s_BlastAaExtendRight, of the hits whose left extension reached the first word.
		rounds = 0;
		for (int l = 0; l < m; l++) {
			leftScore[l] = maxScore[l];
			leftD[l] = n[l] - best[l] + 1;
			right[l] = leftD[l] >= sR[l] - sLeft[l];
			rightScore[l] = 0;
			rightD[l] = 0;
			last[l] = sR[l];
			if (right[l]) {
				n[l] = Math.min(length - sR[l], finder[l].queryLength() - qR[l]);
				best[l] = -1;
				score[l] = maxScore[l];
				floor[l] = Math.max(0, maxScore[l] - finder[l].xDrop());
				i[l] = 0;
				rounds = Math.max(rounds, n[l]);
			} else {
				n[l] = 0;
				i[l] = 0;
			}
		}
		for (int k = 0; k < rounds; k++) {
			boolean any = false;
			for (int l = 0; l < m; l++) {
				final int at = i[l];
				if (at < n[l]) {
					final int sc = score[l] + WordFinder.PADDED[(finder[l].rows()[qR[l] + at]
							+ residues[start + sR[l] + at]) & 1023];
					score[l] = sc;
					if (sc > maxScore[l]) {
						maxScore[l] = sc;
						best[l] = at;
						floor[l] = Math.max(0, sc - finder[l].xDrop());
					}
					// At a break the extension's end is this step; n marks the lane done.
					if (sc <= floor[l]) {
						n[l] = at;
					} else {
						i[l] = at + 1;
					}
					any = true;
				}
			}
			if (!any) {
				break;
			}
		}
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
