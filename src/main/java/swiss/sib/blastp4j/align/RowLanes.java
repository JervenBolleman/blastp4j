package swiss.sib.blastp4j.align;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * {@link GappedAligner#semiGapped}'s rows, several cells of one row at a time. Step 8 of the acceleration plan: on
 * AArch64 the lanes of {@link SemiGappedLanes} are off, so the score-only extensions run one cell at a time, and the
 * rows are wide enough to be worth vectorising (the real extensions of CDK1_HUMAN against Swiss-Prot average 38.5
 * cells of band over 3.08 M rows).
 * <p>
 * The cells are 16 bits wide, as BLAST's own SSE code, SWIPE and parasail all make them: a row of 38.5 cells is then
 * 2.4 groups of 16 lanes rather than 4.8 of 8, which halves both the bookkeeping a row pays and the operations a cell
 * pays. The price is a range check, and {@link #rows} hands an extension that could leave the range to the cells.
 * <p>
 * Three things make a row look unvectorisable, and none of them is:
 * <ul>
 * <li><b>The gap along the row.</b> {@code G(j) = max(G(j - 1) - ext, H(j - 1) - goe)} is a serial chain, but since
 * {@code ext <= goe} a gap routed through an earlier {@code G} pays {@code goe} twice and never wins, so the chain is
 * the prefix maximum {@code A(t) = max(k <= t)(H0(k) + k * ext)} over the scores <em>before</em> the gap
 * ({@code H0 = max(diagonal, gap down the column)}). {@code log2(lanes)} shift-and-maximum steps replace it, and one
 * value carries from one group of lanes to the next. Pinned by {@code RowScanTest}.</li>
 * <li><b>The X-drop.</b> BLAST freezes both gap values at a dropped cell and tightens the floor within the row, which
 * would make every cell depend on the one before it through the branch as well. It does not: a dropped cell's score is
 * below the floor, the floor never falls, and a value below the floor only ever wins a maximum against values below
 * the floor, so both freezes are unobservable and the scores do not depend on the drop decisions at all. The row is
 * therefore computed without the floor and the X-drop is a mask over it, against a floor that is itself a prefix
 * maximum. Pinned by {@code UniformBandTest}.</li>
 * <li><b>The substitution scores.</b> A gather is what made four lanes across extensions lose (fff62be: one
 * {@code ld1w} was 11.07% of the program's samples). A matrix row is {@link Blosum62#SIZE} bytes and BLOSUM62's scores
 * fit a byte, so two table registers hold 32 and one two-vector {@code selectFrom} — {@code TBL} on AArch64, a pair of
 * {@code vpshufb} on x86-64 — reads them all. {@link #of} returns {@link GappedAligner.Rows} for a matrix that does
 * not fit a byte (the composition-based stage scales by 32), so this class may assume it does.</li>
 * </ul>
 * The band's bookkeeping stays scalar, because it is per row and not per cell: the leading run of dropped cells, the
 * last surviving cell, and the short growth loop at the right end.
 */
final class RowLanes implements GappedAligner.RowLoops {
	private static final VectorSpecies<Short> S = ShortVector.SPECIES_PREFERRED;
	/** One 128-bit register per half of the table a {@code selectFrom} reads: 32 bytes, a padded matrix row. */
	private static final VectorSpecies<Byte> B = ByteVector.SPECIES_128;
	private static final int L = S.length();
	/** Lane 0, 1, 2 ... of a group, for the prefix maximum's lifted scores. */
	private static final ShortVector IOTA = ShortVector.zero(S).addIndex(1);
	/**
	 * The no-gap sentinel. Not {@link GappedAligner#MININT}: the prefix-maximum form adds {@code k * ext} to it, which
	 * overflows, and in 16 bits there is no room for a value that large anyway. Low enough that no live cell reaches
	 * it ({@link #dropLimit} keeps {@code xDrop} below it) and high enough to leave the arithmetic room below.
	 */
	private static final short FLOOR = -(1 << 13);
	private static final int N = Blosum62.SIZE;
	/** The sentinels a shift brings in. */
	private static final ShortVector SENTINEL = ShortVector.broadcast(S, FLOOR);
	/**
	 * The lanes moved up by 1, 2, 4 and 8, as shuffles whose low lanes are out of range, so that
	 * {@code v.rearrange(UP, SENTINEL)} takes them from {@link #SENTINEL}. {@code slice}, which is the API's own way to
	 * say this, costs four vector operations where this costs one: it compares an iota against the origin, rearranges
	 * both vectors and blends. The ones a shorter vector has no lanes for are never read; {@link #L} is a static final,
	 * so C2 folds the tests away.
	 */
	private static final VectorShuffle<Short> UP1 = VectorShuffle.iota(S, -1, 1, false),
			UP2 = VectorShuffle.iota(S, -2, 1, false), UP4 = L >= 8 ? VectorShuffle.iota(S, -4, 1, false) : UP1,
			UP8 = L >= 16 ? VectorShuffle.iota(S, -8, 1, false) : UP1;

	private final int goe, ext;
	/** The row loop for the extensions this one hands back, and the two bounds that decide which those are. */
	private final GappedAligner.RowLoops cells;
	private final int bandLimit, dropLimit;
	/** The band, de-interleaved from {@code dp}; see {@link #rows}. Grown together, so one length serves. */
	private short[] rowH = new short[256], rowE = new short[256], rowNext = new short[256];
	/**
	 * What one group of lanes leaves for the next and for the row's bookkeeping: fields rather than returned values, so
	 * that no vector crosses a method boundary. One instance per aligner, so one per thread, as {@link GappedAligner}
	 * is.
	 * <p>
	 * <b>What decides whether C2 emits this kernel as vectors or calls it.</b> Not the number of operations in one
	 * method: 123 of them in a row still compile to vectors (measured). C2 intrinsifies a Vector API operation only
	 * when the vector class and the lane count it passes to {@code VectorSupport} are compile-time constants, and only
	 * within one compilation; where that fails the operation runs its Java implementation, a loop over an array, and
	 * allocates. Three things broke it here, each worth 100 to 300 ms an invocation of {@code RowKernelBench}:
	 * <ul>
	 * <li>a <b>variable lane index</b> ({@code lane(lanes - 1)} for the last, partial group of a row). C2 then needs
	 * the vector as a real object on that path, and with it every vector feeding it, so the whole group fell back:
	 * 9.5 GB of {@code int[]} an invocation. The lanes past the band carry {@link #FLOOR} instead, so every lane is
	 * read at the constant {@code L - 1};</li>
	 * <li>a <b>vector passed as an argument</b>: only a receiver is profiled, so the table registers arrived as the
	 * abstract declared type and the table lookup — the point of the kernel — had no constant lane count. They are
	 * loaded inside {@link #group} instead;</li>
	 * <li>a <b>shared helper called twice</b>: C2 inlined the prefix maximum at the first call site and refused it at
	 * the second, and the vector crossing that boundary cost 7 of 29 young collections. Both are written out in
	 * {@link #group}.</li>
	 * </ul>
	 */
	private int gapRow, bestScore, floor, last, newFirst;
	private boolean leading;

	private RowLanes(int goe, int ext, GappedAligner.RowLoops cells, int maxCell) {
		this.goe = goe;
		this.ext = ext;
		this.cells = cells;
		// An alignment through a cell of row i and column j has at most min(i, j) aligned pairs and pays for every gap,
		// so no score passes maxCell * min(m, n); the lift of the prefix maximum adds (L - 1) * ext on top of it and a
		// diagonal one more cell.
		this.bandLimit = (Short.MAX_VALUE - (L - 1) * ext - maxCell) / maxCell;
		// A dropped cell keeps the sentinel, and the next row reads it as a diagonal, FLOOR + maxCell. For that to drop
		// in its turn — which is what makes the sentinel behave as BLAST's MININT — it has to be below the lowest floor
		// a row can have, and that is -xDrop, since the best score starts at zero and only grows.
		this.dropLimit = -FLOOR - maxCell - 1;
	}

	/**
	 * The row loop for these gap costs and this matrix: this class where the vectors can be used, and null where they
	 * cannot, so that {@link GappedAligner} keeps {@code cells} ({@link GappedAligner.Rows}). {@code matrix} is the
	 * aligner's, {@link Blosum62#SIZE} wide; {@code cells} is where an extension too long or too wide for 16 bits goes.
	 */
	static GappedAligner.RowLoops of(int goe, int ext, int[] matrix, GappedAligner.RowLoops cells) {
		if (L < 4 || ext > goe || !fits(matrix)) {
			return null;
		}
		int maxCell = 1;
		for (int score : matrix) {
			maxCell = Math.max(maxCell, score);
		}
		// Room below the sentinel for everything a cell at the sentinel can reach. Nothing drifts below it row by row:
		// a cell is either dropped, and then it holds the sentinel itself, or live, and then it is at least the floor,
		// so the lowest value anywhere is a gap opened from a sentinel and carried across one group of lanes. Four gap
		// openings and two groups of extension is generous for that, and these costs never come near it.
		return FLOOR - 4 * goe - 2 * L * ext - maxCell > Short.MIN_VALUE
				? new RowLanes(goe, ext, cells, maxCell)
				: null;
	}

	/** Whether this row loop can be used for {@code matrix}: its scores have to fit a byte, as the lookup is one. */
	static boolean fits(int[] matrix) {
		for (int score : matrix) {
			if (score != (byte) score) {
				return false;
			}
		}
		return true;
	}

	@Override
	public int rows(GappedAligner g, byte[] a, int aBase, int m, int n, boolean reverse, int xDrop, int bSize) {
		// What 16 bits cannot hold goes to the cells. Both bounds are per extension, not per cell: a score cannot pass
		// maxCell * min(m, n), and the sentinel is only BLAST's MININT while it stays below every floor.
		if (xDrop > dropLimit || Math.min(m, n) > bandLimit) {
			return cells.rows(g, a, aBase, m, n, reverse, xDrop, bSize);
		}
		final int goe = this.goe, ext = this.ext;
		final int[] dp = g.dp, scores = g.matrix.scores;
		if (rowH.length < n + 2 + GappedAligner.PAD) {
			final int size = Math.max(n + 2 + GappedAligner.PAD, rowH.length * 2);
			rowH = new short[size];
			rowE = new short[size];
			rowNext = new short[size];
		}
		// The band, de-interleaved from dp once per extension: a vector of scores cannot be read from dp's interleaved
		// layout at all. Column j is h[j + 1] and e[j + 1], so that the diagonal's load from column j - 1 is in bounds
		// at the band's first column; the arrays hold n + 2 + PAD columns, and the highest a row reaches is n + 1.
		//
		// The scores go in two arrays, read and written in turn, because the groups write where the diagonal reads,
		// one column apart: writing in place would destroy the row above's score that the next group needs, and a
		// second array is cheaper than carrying that value, which costs a lane insert in every group. The gap down a
		// column is read and written at the same column, so it stays in one. Columns outside the band hold what a row
		// two back left there and are never read: the band's first column takes the sentinel for its diagonal, and
		// the growth below writes every column past its last one.
		short[] h = rowH, e = rowE, next = rowNext;
		final byte[] bs = g.bLetters;
		for (int j = 0; j < bSize; j++) {
			h[j + 1] = (short) dp[2 * j];
			e[j + 1] = (short) dp[2 * j + 1];
		}
		// A matrix row as two table registers, padded to 32 bytes.
		final byte[] table = new byte[32];
		bestScore = 0;
		floor = -xDrop;
		int first = 0;
		// The diagonal of the band's first column, which an earlier row dropped: the sentinel at h[first] (column
		// first - 1) rather than a blend into lane 0 of the first group, so that the group has no vector of its own
		// that a branch decides. The group loop writes next[c + 1] onwards, so it never overwrites it.
		h[0] = FLOOR;
		for (int ai = 1; ai <= m; ai++) {
			final int letter = GappedAligner.letterA(a, aBase, m, ai, reverse);
			for (int k = 0; k < N; k++) {
				table[k] = (byte) scores[letter * N + k];
			}
			gapRow = FLOOR;
			last = first;
			newFirst = first;
			leading = true;
			for (int c = first; c < bSize; c += L) {
				group(g, h, next, e, bs, table, c, bSize, ai, xDrop);
			}
			final short[] written = next;
			next = h;
			h = written;
			first = newFirst;
			if (first == bSize) {
				break;
			}
			h[first] = FLOOR;
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= floor && bSize <= n) {
					h[bSize + 1] = (short) gapRow;
					e[bSize + 1] = (short) (gapRow - goe);
					gapRow -= ext;
					bSize++;
				}
			}
			if (bSize <= n) {
				h[bSize + 1] = FLOOR;
				e[bSize + 1] = FLOOR;
				bSize++;
			}
		}
		return bestScore;
	}

	/**
	 * One group of lanes of a row: the scores of its cells, the gap along the row through them, and the X-drop as a
	 * mask over them. Reads and writes the fields above, which is what keeps every vector inside one compilation.
	 */
	private void group(GappedAligner g, short[] h, short[] next, short[] e, byte[] bs, byte[] table, int c, int bSize,
			int ai, int xDrop) {
		final short goe = (short) this.goe, ext = (short) this.ext;
		// The two table registers are loaded here and not once per row: passed in, they are parameters of the declared
		// type ByteVector, C2 has no exact type for them (only a receiver is profiled, not an argument), so the vector
		// class and the lane count that selectFromTwoVectorOp needs as constants are not constants, and the table
		// lookup — the whole point of this kernel — runs in Java and allocates. Two L1 loads a group is the price of
		// the lookup being a TBL.
		final ByteVector low = ByteVector.fromArray(B, table, 0), high = ByteVector.fromArray(B, table, 16);
		final ShortVector step = IOTA.mul(ext);
		final int lanes = Math.min(L, bSize - c);
		// The lanes past the band's last column. The loads and the stores need no mask of their own: the arrays hold
		// GappedAligner.PAD columns beyond the band, so a lane past it reads what an earlier row left and writes where
		// nothing reads before the band grows there.
		final VectorMask<Short> out = IOTA.compare(VectorOperators.GE, (short) lanes);
		// H0 = max(diagonal, gap down the column). B's letter of the diagonal is the one before the column, as the row
		// loop reads it, and its score comes from the table rather than from a gather.
		final int bOff = Math.max(c - 1, 0);
		final ByteVector letters = ByteVector.fromArray(B, bs, bOff).and((byte) 31);
		ShortVector subst = (ShortVector) letters.selectFrom(low, high).convertShape(VectorOperators.B2S, S, 0);
		if (c == 0) {
			// The letters start at the band's column 0, one later than the lanes want; lane 0's diagonal is the
			// sentinel the row loop left at h[first] in any case.
			subst = subst.rearrange(UP1, SENTINEL);
		}
		final ShortVector e0 = ShortVector.fromArray(S, e, c + 1);
		// Past the band the sentinel, so that the two prefix maximums below carry the band's own last lane into lane
		// L - 1 and every lane this group reads is read at a constant index.
		final ShortVector h0 = ShortVector.fromArray(S, h, c).add(subst).max(e0).blend(FLOOR, out);
		// The gap along the row as a prefix maximum: A(t) = max(k <= t)(H0(k) + k * ext), then
		// G(t) = max(carried - t * ext, A(t - 1) - goe - (t - 1) * ext) and H = max(H0, G).
		ShortVector lifted = h0.add(step);
		lifted = lifted.max(lifted.rearrange(UP1, SENTINEL));
		lifted = lifted.max(lifted.rearrange(UP2, SENTINEL));
		if (L >= 8) {
			lifted = lifted.max(lifted.rearrange(UP4, SENTINEL));
		}
		if (L >= 16) {
			lifted = lifted.max(lifted.rearrange(UP8, SENTINEL));
		}
		final ShortVector gap = ShortVector.broadcast(S, (short) gapRow).sub(step)
				.max(lifted.rearrange(UP1, SENTINEL).sub(step).add((short) (ext - goe)));
		final ShortVector hv = h0.max(gap);
		// What the next group carries, from the lifted prefix maximum rather than from the last cell:
		// G(c + lanes) = max(G(c) - lanes * ext, A - goe - (lanes - 1) * ext).
		final int a = lifted.lane(L - 1);
		gapRow = Math.max(gapRow - lanes * ext, a - goe - (lanes - 1) * ext);
		// The X-drop, as a mask against a floor that is a prefix maximum of the scores. Including a dropped cell in
		// that maximum is harmless: it is below the floor, so it is never the maximum. The lanes past the band are at
		// the sentinel, so they drop with it and need no mask of their own.
		ShortVector run = hv.max(hv.rearrange(UP1, SENTINEL));
		run = run.max(run.rearrange(UP2, SENTINEL));
		if (L >= 8) {
			run = run.max(run.rearrange(UP4, SENTINEL));
		}
		if (L >= 16) {
			run = run.max(run.rearrange(UP8, SENTINEL));
		}
		final VectorMask<Short> dropped = hv.lt(run.max((short) bestScore).sub((short) xDrop));
		// A dropped cell's score is killed so that it cannot feed a diagonal into the next row.
		hv.blend(FLOOR, dropped).intoArray(next, c + 1);
		hv.sub(goe).max(e0.sub(ext)).intoArray(e, c + 1);
		// The bookkeeping, from the mask's bits: one read of it rather than three.
		final VectorMask<Short> alive = dropped.not();
		final long live = alive.toLong();
		if (live != 0) {
			last = c + 63 - Long.numberOfLeadingZeros(live);
			// The group's best, read from the prefix maximum's last lane. A group all of whose cells are dropped
			// cannot beat the best score (they are all below the floor), so the lanes need no mask here either.
			final int top = run.lane(L - 1);
			if (top > bestScore) {
				bestScore = top;
				floor = top - xDrop;
				g.aOffset = ai;
				g.bOffset = c + hv.eq((short) top).and(alive).firstTrue();
			}
			if (leading) {
				newFirst = c + Long.numberOfTrailingZeros(live);
				leading = false;
			}
		} else if (leading) {
			newFirst = c + lanes;
		}
	}
}
