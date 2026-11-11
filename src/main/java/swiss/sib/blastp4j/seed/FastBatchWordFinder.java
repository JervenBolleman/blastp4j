package swiss.sib.blastp4j.seed;

import java.util.List;

/**
 * The word finder of {@code blastp -task blastp-fast} for a batch of queries: {@link BatchWordFinder} with BLAST's
 * compressed-alphabet scan ({@code s_BlastCompressedAaScanSubject}) and its {@link CompressedTable}. A subject position
 * whose 5 letters are all in the {@link CompressedAlphabet} is looked up; the others (with {@code X}, {@code U},
 * {@code *}, {@code O} or {@code -}) are skipped, as BLAST skips them. The two-hit logic is {@link BatchWordFinder}'s
 * with words of 5 letters, and the hits that are not just stored go to the query's {@link FastWordFinder}.
 * <p>
 * Not thread safe: use one per thread; the {@link CompressedTable} can be shared.
 */
public final class FastBatchWordFinder implements BatchScan {
	private static final int WORD = CompressedTable.WORD;
	/** The weight of a word's newest letter: 15^4. */
	private static final int NEWEST = 50625;
	/** How many positions with hits ahead the scan touches the cells of. */
	private static final int AHEAD = 16;

	private final CompressedTable table;
	private final FastWordFinder[] finders;
	/** The diagonals' states in 16 bits, as {@link BatchWordFinder}'s. */
	private final char[] diags;
	/** The diagonals as ints, for a subject too long for 16 bits. */
	private int[] wide;
	private static final int LIMIT = 1 << 15, EXTENDED = 1 << 15;
	private final int window, diagMask;
	private int diagOffset;
	/** The positions with hits of the subject being scanned, and their words. */
	private int[] words = new int[0], at = new int[0];
	/** The sum of the cells touched ahead, kept so that the loads stay. */
	private int sink;

	public FastBatchWordFinder(byte[][] queries, SeedParameters[] parameters, CompressedTable table) {
		this.table = table;
		window = parameters[0].window();
		diags = new char[table.size()];
		diagMask = table.size() - 1;
		finders = new FastWordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new FastWordFinder(queries[i], parameters[i]);
		}
		diagOffset = window;
	}

	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		if (length + window >= LIMIT) {
			findLong(residues, start, length, out);
			return;
		}
		if (diagOffset + length >= LIMIT) {
			java.util.Arrays.fill(diags, (char) 0);
			diagOffset = window;
		}
		final int[] cells = table.cells(), entry = table.entry(), query = table.query(), base = table.base();
		final long[] present = table.present();
		final char[] d = diags;
		final int diagMask = this.diagMask & (d.length - 1);
		final int offset = diagOffset, windowL = window, lastL = window - WORD - 1;
		final byte[] compress = CompressedAlphabet.COMPRESS;
		final FastWordFinder[] f = finders;
		// First the positions whose cell has hits, with their words; then their hits, each position touching the cells
		// of the one AHEAD further on, so that the misses on cells (6 MB) overlap rather than wait one after the other.
		if (words.length < length) {
			words = new int[length];
			at = new int[length];
		}
		final int[] w = words, pos = at;
		int n = 0;
		// The word ending at each position, its oldest letter lowest; good counts the letters since one outside the
		// alphabet.
		int word = 0, good = 0;
		for (int p = 0; p < length; p++) {
			final int c = compress[residues[start + p]];
			if (c >= CompressedAlphabet.SIZE) {
				good = 0;
				word = 0;
				continue;
			}
			word = word / CompressedAlphabet.SIZE + c * NEWEST;
			if (++good < WORD) {
				continue;
			}
			if ((present[word >>> 6] >>> word & 1) == 0) {
				continue; // no hits: cells not read
			}
			w[n] = word;
			pos[n++] = p - (WORD - 1);
		}
		int touched = 0;
		for (int j = 0; j < n; j++) {
			if (j + AHEAD < n) {
				touched += cells[2 * w[j + AHEAD]];
			}
			final int cell = 2 * w[j], s = pos[j];
			// The cell's first entry from its own line of cells, the others from entry. A counted loop, which C2
			// keeps without a safepoint poll per hit: each hit loads the next one's entry (entry has one more).
			int e = cells[cell + 1];
			for (int k = cells[cell], end = cells[cell + 2]; k < end; k++) {
				// As BatchWordFinder.find.
				final int diag = (e - s) & diagMask;
				final int state = d[diag];
				final int diff = s - (state - offset), x = diff - WORD;
				if ((~(x | (lastL - x)) | state << 16) >= 0) {
					final int older = (diff - windowL) >> 31;
					d[diag] = (char) ((state & older) | ((s + offset) & ~older));
				} else {
					final int i = query[k];
					final int next = f[i].next(residues, start, length, e - base[i], s, out[i],
							(state & (EXTENDED - 1)) | (state & EXTENDED) << 15, offset);
					// The diagonal again from the entry (entry[k] is the cell's first one too), so that it need not
					// be kept across the call.
					d[(entry[k] - s) & diagMask] = (char) (next | next >>> 15);
				}
				e = entry[k + 1];
			}
		}
		sink += touched;
		diagOffset += length + window;
	}

	/**
	 * {@link #find} for a subject whose positions do not fit 16 bits, with diagonals of ints, cleared for it. The
	 * 16-bit states and their offset stay as they are.
	 */
	void findLong(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cells = table.cells(), entry = table.entry(), query = table.query(), base = table.base();
		final long[] present = table.present();
		final byte[] compress = CompressedAlphabet.COMPRESS;
		if (wide == null) {
			wide = new int[diags.length];
		} else {
			java.util.Arrays.fill(wide, 0);
		}
		final int[] d = wide;
		final int offset = window;
		int word = 0, good = 0;
		for (int p = 0; p < length; p++) {
			final int c = compress[residues[start + p]];
			if (c >= CompressedAlphabet.SIZE) {
				good = 0;
				word = 0;
				continue;
			}
			word = word / CompressedAlphabet.SIZE + c * NEWEST;
			if (++good < WORD || (present[word >>> 6] >>> word & 1) == 0) {
				continue;
			}
			final int cell = 2 * word, s = p - (WORD - 1);
			int e = cells[cell + 1];
			for (int k = cells[cell], end = cells[cell + 2]; k < end; k++) {
				final int i = query[k], diag = (e - s) & (d.length - 1);
				d[diag] = finders[i].next(residues, start, length, e - base[i], s, out[i], d[diag], offset);
				e = entry[k + 1];
			}
		}
	}
}
