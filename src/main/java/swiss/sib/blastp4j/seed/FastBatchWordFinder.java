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

	private final CompressedTable table;
	private final FastWordFinder[] finders;
	private final int[] diags;
	private final int window, diagMask;
	private int diagOffset;

	public FastBatchWordFinder(byte[][] queries, SeedParameters[] parameters, CompressedTable table) {
		this.table = table;
		window = parameters[0].window();
		diags = new int[table.size()];
		diagMask = table.size() - 1;
		finders = new FastWordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new FastWordFinder(queries[i], parameters[i], diags, table.base()[i]);
		}
		diagOffset = window;
	}

	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cells = table.cells(), entry = table.entry(), query = table.query(), base = table.base();
		final long[] present = table.present();
		final int[] d = diags;
		final int diagMask = this.diagMask & (d.length - 1);
		final int offset = diagOffset, windowL = window, lastL = window - WORD - 1;
		final byte[] compress = CompressedAlphabet.COMPRESS;
		final FastWordFinder[] f = finders;
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
				continue; // no hits: cells, 6 MB, not read
			}
			final int s = p - (WORD - 1);
			// The cell's first entry from its own line of cells, the others from entry.
			int k = cells[2 * word], e = cells[2 * word + 1];
			for (final int end = cells[2 * word + 2];; e = entry[k]) {
				// As BatchWordFinder.find.
				final int diag = (e - s) & diagMask;
				final int state = d[diag];
				final int diff = s - (state - offset), x = diff - WORD;
				if ((~(x | (lastL - x)) | state << 1) >= 0) {
					final int older = (diff - windowL) >> 31;
					d[diag] = (state & older) | ((s + offset) & ~older);
				} else {
					final int i = query[k];
					f[i].hit(residues, start, length, e - base[i], s, out[i]);
				}
				if (++k >= end) {
					break;
				}
			}
		}
		for (FastWordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}
}
