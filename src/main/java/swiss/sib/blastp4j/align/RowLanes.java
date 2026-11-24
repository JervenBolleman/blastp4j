package swiss.sib.blastp4j.align;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * {@link GappedAligner#semiGapped}'s rows, several cells of one row at a time. Step 8 of the acceleration plan: on
 * AArch64 the lanes of {@link SemiGappedLanes} are off, so the score-only extensions run one cell at a time, and the
 * rows are wide enough to be worth vectorising (the real extensions of CDK1_HUMAN against Swiss-Prot average 38.5
 * cells of band over 3.08 M rows).
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
	private static final VectorSpecies<Integer> I = IntVector.SPECIES_PREFERRED;
	/** One 128-bit register per half of the table a {@code selectFrom} reads: 32 bytes, a padded matrix row. */
	private static final VectorSpecies<Byte> B = ByteVector.SPECIES_128;
	private static final int L = I.length();
	/** Lane 0, 1, 2 ... of a group, for the prefix maximum's lifted scores. */
	private static final IntVector IOTA = IntVector.zero(I).addIndex(1);
	/**
	 * The no-gap sentinel. Not {@link GappedAligner#MININT}: the prefix-maximum form adds {@code k * ext} to it, which
	 * overflows. Far enough below any reachable score that it never wins a maximum.
	 */
	private static final int FLOOR = -(1 << 28);
	private static final int N = Blosum62.SIZE;
	/** Lane 0 alone, and a vector of sentinels to shift in. */
	private static final VectorMask<Integer> LANE0 = VectorMask.fromLong(I, 1);
	private static final IntVector SENTINEL = IntVector.broadcast(I, FLOOR);

	private final int goe, ext;
	/**
	 * What one group of lanes leaves for the next and for the row's bookkeeping. Fields, not returned values, because
	 * {@link #group} has to stay small enough for C2 to inline the vector helpers into it: at 932 bytes of bytecode in
	 * one method they hit {@code NodeCountInliningCutoff}, {@code slice} fell back to its Java implementation, and the
	 * kernel allocated a vector per call (3.5 times slower than the cells, 55 young collections against 20). One
	 * instance per aligner, so one per thread, as {@link GappedAligner} is.
	 */
	private int gapRow, bestScore, floor, last, newFirst;
	private boolean leading;

	private RowLanes(int goe, int ext) {
		this.goe = goe;
		this.ext = ext;
	}

	/**
	 * The row loop for these gap costs and this matrix: this class where the vectors can be used, and null where they
	 * cannot, so that {@link GappedAligner} keeps {@link GappedAligner.Rows}. {@code matrix} is the aligner's,
	 * {@link Blosum62#SIZE} wide.
	 */
	static GappedAligner.RowLoops of(int goe, int ext, int[] matrix) {
		return L >= 4 && ext <= goe && fits(matrix) ? new RowLanes(goe, ext) : null;
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
		final int goe = this.goe, ext = this.ext;
		final int[] dp = g.dp, scores = g.matrix.scores;
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
		int[] h = g.best, e = g.bestGap, next = g.rowNext;
		final byte[] bs = g.bLetters;
		for (int j = 0; j < bSize; j++) {
			h[j + 1] = dp[2 * j];
			e[j + 1] = dp[2 * j + 1];
		}
		// A matrix row as two table registers, padded to 32 bytes.
		final byte[] table = new byte[32];
		bestScore = 0;
		floor = -xDrop;
		int first = 0;
		for (int ai = 1; ai <= m; ai++) {
			final int letter = GappedAligner.letterA(a, aBase, m, ai, reverse);
			for (int k = 0; k < N; k++) {
				table[k] = (byte) scores[letter * N + k];
			}
			final ByteVector low = ByteVector.fromArray(B, table, 0), high = ByteVector.fromArray(B, table, 16);
			gapRow = FLOOR;
			last = first;
			newFirst = first;
			leading = true;
			for (int c = first; c < bSize; c += L) {
				group(g, h, next, e, bs, low, high, c, bSize, first, ai, xDrop);
			}
			final int[] written = next;
			next = h;
			h = written;
			first = newFirst;
			if (first == bSize) {
				break;
			}
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= floor && bSize <= n) {
					h[bSize + 1] = gapRow;
					e[bSize + 1] = gapRow - goe;
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
	 * mask over them. Reads and writes the fields above, which is what keeps it small.
	 */
	private void group(GappedAligner g, int[] h, int[] next, int[] e, byte[] bs, ByteVector low, ByteVector high, int c,
			int bSize, int first, int ai, int xDrop) {
		final int goe = this.goe, ext = this.ext;
		final IntVector step = IOTA.mul(ext), stepBefore = IOTA.sub(1).mul(ext);
		final int lanes = Math.min(L, bSize - c);
		// Which lanes are in the band. The loads and the stores need no mask of their own: the arrays hold
		// GappedAligner.PAD columns beyond the band, a lane past its last column reads what an earlier row left and
		// writes where nothing reads, and the prefix maximums below only ever move a lane up, so those lanes cannot
		// reach the ones in the band.
		final VectorMask<Integer> in = IOTA.lt(lanes);
		// H0 = max(diagonal, gap down the column). B's letter of the diagonal is the one before the column, as the row
		// loop reads it, and its score comes from the table rather than from a gather. The band's first column has no
		// diagonal: an earlier row dropped the cell above-left of it.
		final int bOff = Math.max(c - 1, 0);
		final ByteVector letters = ByteVector.fromArray(B, bs, bOff).and((byte) 31);
		IntVector subst = (IntVector) letters.selectFrom(low, high).convertShape(VectorOperators.B2I, I, 0);
		if (c == 0) {
			// The letters start at the band's column 0, one later than the lanes want; lane 0's diagonal is the
			// sentinel below in any case.
			subst = shiftOne(subst);
		}
		IntVector diag = IntVector.fromArray(I, h, c).add(subst);
		if (c == first) {
			diag = diag.blend(FLOOR, LANE0);
		}
		final IntVector e0 = IntVector.fromArray(I, e, c + 1);
		final IntVector h0 = diag.max(e0);
		// The gap along the row as a prefix maximum: A(t) = max(k <= t)(H0(k) + k * ext), then
		// G(t) = max(carried - t * ext, A(t - 1) - goe - (t - 1) * ext) and H = max(H0, G).
		final IntVector lifted = prefixMax(h0.add(step));
		final IntVector gap = IntVector.broadcast(I, gapRow).sub(step)
				.max(shiftOne(lifted).sub(stepBefore).sub(goe));
		final IntVector hv = h0.max(gap);
		// What the next group carries, from the lifted prefix maximum rather than from the last cell, so that the lane
		// is a constant wherever the group is full: G(c + lanes) = max(G(c) - lanes * ext, A - goe - (lanes - 1) * ext).
		final int a = lanes == L ? lifted.lane(L - 1) : lifted.lane(lanes - 1);
		gapRow = Math.max(gapRow - lanes * ext, a - goe - (lanes - 1) * ext);
		// The X-drop, as a mask against a floor that is a prefix maximum of the scores. Including a dropped cell in
		// that maximum is harmless: it is below the floor, so it is never the maximum.
		final IntVector run = prefixMax(hv);
		final VectorMask<Integer> dropped = hv.lt(run.max(bestScore).sub(xDrop)).and(in);
		// A dropped cell's score is killed so that it cannot feed a diagonal into the next row.
		hv.blend(FLOOR, dropped).intoArray(next, c + 1);
		hv.sub(goe).max(e0.sub(ext)).intoArray(e, c + 1);
		// The bookkeeping, from the mask's bits: one read of it rather than three.
		final long live = dropped.not().and(in).toLong();
		if (live != 0) {
			last = c + 63 - Long.numberOfLeadingZeros(live);
			// The group's best, read from the prefix maximum's last lane. A group all of whose cells are dropped
			// cannot beat the best score (they are all below the floor), so the lanes need no mask here.
			final int top = lanes == L ? run.lane(L - 1) : run.lane(lanes - 1);
			if (top > bestScore) {
				bestScore = top;
				floor = top - xDrop;
				g.aOffset = ai;
				g.bOffset = c + hv.eq(top).and(VectorMask.fromLong(I, live)).firstTrue();
			}
			if (leading) {
				newFirst = c + Long.numberOfTrailingZeros(live);
				leading = false;
			}
		} else if (leading) {
			newFirst = c + lanes;
		}
	}

	/**
	 * {@code v} moved up by one lane, the lane below it filled with {@link #FLOOR}.
	 * <p>
	 * {@code slice}'s origin has to be a constant: with a loop's variable there, C2 builds a shuffle object at run time
	 * instead of emitting the shift, which allocates and was most of this kernel's first cost. That is why the prefix
	 * maximum below is unrolled and why {@link #L} is a static final the folds depend on.
	 */
	private static IntVector shiftOne(IntVector v) {
		return SENTINEL.slice(L - 1, v);
	}

	/** The inclusive prefix maximum over the lanes: {@code log2(L)} shift-and-maximum steps, each shift a constant. */
	private static IntVector prefixMax(IntVector v) {
		v = v.max(SENTINEL.slice(L - 1, v));
		v = v.max(SENTINEL.slice(L - 2, v));
		if (L >= 8) {
			v = v.max(SENTINEL.slice(L - 4, v));
		}
		if (L >= 16) {
			v = v.max(SENTINEL.slice(L - 8, v));
		}
		return v;
	}
}
