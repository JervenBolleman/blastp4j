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
	 * The queries' lookup tables merged: the hits of word {@code w} are entries {@code [start[w], start[w + 1])}, query
	 * {@code query[k]} at offset {@code offset[k]}, by query and then as in that query's table.
	 */
	public record Table(int[] start, int[] query, int[] offset) {
		public static Table of(LookupTable[] tables) {
			int total = 0;
			for (LookupTable t : tables) {
				total += t.offsets().length;
			}
			final int[] start = new int[LookupTable.WORDS + 1], query = new int[total], offset = new int[total];
			int k = 0;
			for (int w = 0; w < LookupTable.WORDS; w++) {
				start[w] = k;
				for (int i = 0; i < tables.length; i++) {
					final int[] ts = tables[i].start(), to = tables[i].offsets();
					for (int j = ts[w]; j < ts[w + 1]; j++) {
						query[k] = i;
						offset[k] = to[j];
						k++;
					}
				}
			}
			start[LookupTable.WORDS] = k;
			return new Table(start, query, offset);
		}
	}

	private final Table table;
	private final WordFinder[] finders;

	/** The queries, their lookup tables and seeding parameters, by query; {@code table} merges the tables. */
	public BatchWordFinder(byte[][] queries, LookupTable[] tables, SeedParameters[] parameters, Table table) {
		this.table = table;
		finders = new WordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			finders[i] = new WordFinder(queries[i], tables[i], parameters[i]);
		}
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} for each query, adding
	 * those of query {@code i} to {@code out[i]}, as {@link WordFinder#find} does for one query.
	 */
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cellStart = table.start(), query = table.query(), offset = table.offset();
		final WordFinder[] f = finders;
		final int mask = LookupTable.WORDS - 1;
		if (length >= LookupTable.WORD_SIZE) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - LookupTable.WORD_SIZE; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					final int i = query[k];
					f[i].probe(residues, start, length, offset[k], s, out[i]);
				}
			}
		}
		for (WordFinder finder : f) {
			finder.endSubject(length);
		}
	}
}
