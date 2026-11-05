package swiss.sib.blastp4j.seed;

import java.util.List;

/**
 * The word finder for a batch of queries: BLAST scans each subject once for a batch of queries rather than once per
 * query. Each subject position's word is looked up once, in a table that merges the queries' lookup tables, and each
 * hit goes to its query's {@link WordFinder}, with that query's diagonal array, X-drop and cutoff. A query's hits
 * arrive in the order its own word finder takes them (by subject position, then by query offset), so each query gets
 * the alignments {@link WordFinder#find} gives it alone.
 * <p>
 * Not thread safe: use one per thread; the {@link Table} can be shared.
 */
public final class BatchWordFinder {
	/**
	 * The queries' lookup tables merged, with where each query's diagonals are in the batch's diagonal array: the hits
	 * of word {@code w} are entries {@code [start[w], start[w + 1])}, by query and then as in that query's table.
	 * <p>
	 * Each query has a region of {@code 2^b} diagonals ({@link WordFinder#diagonals}), starting at a multiple of its
	 * size. An entry is {@code b << POSITION_BITS | (region start + query offset)}, so the diagonal of a hit at subject
	 * position {@code s} follows from the entry alone; {@code query[k]} is the query of entry {@code k}, read only for
	 * hits that go on to the two-hit logic. {@code base[i]} is the start of query {@code i}'s region, {@code size} the
	 * length of the array.
	 */
	public record Table(int[] start, int[] entry, int[] query, int[] base, int size) {
		public static Table of(byte[][] queries, LookupTable[] tables, int window) {
			final int n = queries.length;
			// The regions by size, largest first, so that each starts at a multiple of its size.
			final Integer[] order = new Integer[n];
			final int[] bits = new int[n];
			for (int i = 0; i < n; i++) {
				order[i] = i;
				bits[i] = Integer.numberOfTrailingZeros(WordFinder.diagonals(queries[i].length, window));
			}
			java.util.Arrays.sort(order, (a, b) -> Integer.compare(bits[b], bits[a]));
			final int[] base = new int[n];
			long size = 0;
			for (int i : order) {
				base[i] = (int) size;
				size += 1L << bits[i];
			}
			if (size > POSITION_MASK + 1L) {
				throw new IllegalArgumentException("the diagonals of a batch of " + n + " queries (" + size
						+ ") do not fit an entry");
			}
			int total = 0;
			for (LookupTable t : tables) {
				total += t.offsets().length;
			}
			final int[] start = new int[LookupTable.WORDS + 1], entry = new int[total], query = new int[total];
			int k = 0;
			for (int w = 0; w < LookupTable.WORDS; w++) {
				start[w] = k;
				for (int i = 0; i < n; i++) {
					final int[] ts = tables[i].start(), to = tables[i].offsets();
					for (int j = ts[w]; j < ts[w + 1]; j++) {
						query[k] = i;
						entry[k++] = bits[i] << POSITION_BITS | (base[i] + to[j]);
					}
				}
			}
			start[LookupTable.WORDS] = k;
			return new Table(start, entry, query, base, (int) size);
		}
	}

	/** The low bits of an entry: the diagonal region's start plus the query offset. */
	private static final int POSITION_BITS = 26, POSITION_MASK = (1 << POSITION_BITS) - 1;

	private final Table table;
	private final WordFinder[] finders;
	/** The diagonals of all queries, in the regions of {@link Table#base}. */
	private final int[] diags;
	private final int window;
	/** As each finder's {@code diagOffset}, which all move together. */
	private int diagOffset;

	/** The queries, their lookup tables and seeding parameters, by query; {@code table} merges the tables. */
	public BatchWordFinder(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		this.table = table;
		window = parameters[0].window();
		diags = new int[table.size()];
		finders = new WordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new WordFinder(queries[i], tables[i], parameters[i], diags, table.base()[i]);
		}
		diagOffset = window;
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} for each query, adding
	 * those of query {@code i} to {@code out[i]}, as {@link WordFinder#find} does for one query.
	 */
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cellStart = table.start(), entry = table.entry(), query = table.query();
		final int[] d = diags;
		final int offset = diagOffset, windowL = window, lastL = window - LookupTable.WORD_SIZE - 1;
		final WordFinder[] f = finders;
		final int mask = LookupTable.WORDS - 1;
		if (length >= LookupTable.WORD_SIZE) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - LookupTable.WORD_SIZE; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					// The usual hits, decided from the entry: the diagonal's last hit is a window or more back, so
					// this one replaces it, or it overlaps this one, which is then ignored. Only hits to extend and
					// hits on an extended diagonal go to the query's finder.
					final int e = entry[k];
					final int m = (1 << (e >>> POSITION_BITS)) - 1, p = e & POSITION_MASK;
					final int diag = (p & ~m) | ((p - s) & m);
					final int state = d[diag];
					// Without branches but one: x = diff - WORD_SIZE is in [0, last] (a hit to extend) when neither
					// x nor last - x is negative, and the extended flag (bit 30; bit 31 is never set) goes to the
					// sign bit.
					final int diff = s - (state - offset), x = diff - LookupTable.WORD_SIZE;
					if ((~(x | (lastL - x)) | state << 1) >= 0) {
						final int older = (diff - windowL) >> 31; // -1 if the hit overlaps the last one: keep that
						d[diag] = (state & older) | ((s + offset) & ~older);
					} else {
						final int i = query[k];
						f[i].hit(residues, start, length, p & m, s, out[i]);
					}
				}
			}
		}
		for (WordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}
}
