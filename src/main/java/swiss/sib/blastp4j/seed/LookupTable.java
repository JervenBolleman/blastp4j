package swiss.sib.blastp4j.seed;

import static java.util.Arrays.copyOfRange;
import swiss.sib.blastp4j.stats.Blosum62;

/**
 * blastp's lookup table: for every word of three letters, the query positions whose word scores at least the threshold
 * against it (the neighbourhood words), plus the query's own words whatever their score ({@code blast_aalookup.c}).
 * <p>
 * A word is indexed by its three NCBIstdaa codes, five bits each, first letter highest, as in BLAST. The table is
 * stored compressed: the query offsets of word {@code w} are {@code offsets[start[w] .. start[w + 1])}.
 * <p>
 * Instances are immutable.
 */
public record LookupTable(int[] start, int[] offsets, int[] presence) {
	/** Letters per word. */
	public static final int WORD_SIZE = 3;
	/** Bits per letter. */
	static final int CHAR_BITS = 5;
	/** Number of words (some unused, the alphabet has 28 letters). */
	static final int WORDS = 1 << (WORD_SIZE * CHAR_BITS);

	/** Index of the word at {@code pos} of {@code seq}. */
	static int word(byte[] seq, int pos) {
		return (seq[pos] << (2 * CHAR_BITS)) | (seq[pos + 1] << CHAR_BITS) | seq[pos + 2];
	}

	/** The table of a query (NCBIstdaa codes) for neighbourhood score {@code threshold} (11 for blastp). */
	public static LookupTable of(byte[] query, int threshold) {
		final int n = Blosum62.SIZE;
		int[] rowMax = new int[n];
		for (int a = 0; a < n; a++) {
			rowMax[a] = Integer.MIN_VALUE;
			for (int b = 0; b < n; b++) {
				rowMax[a] = Math.max(rowMax[a], Blosum62.score(a, b));
			}
		}
		// Group the query offsets by their word, then find the neighbours once per distinct word.
		int words = Math.max(0, query.length - WORD_SIZE + 1);
		int[] count = new int[WORDS + 1];
		for (int q = 0; q < words; q++) {
			count[word(query, q) + 1]++;
		}
		int[] groupStart = new int[WORDS + 1];
		for (int w = 0; w < WORDS; w++) {
			groupStart[w + 1] = groupStart[w] + count[w + 1];
		}
		int[] grouped = new int[words];
		int[] fill = groupStart.clone();
		for (int q = 0; q < words; q++) {
			grouped[fill[word(query, q)]++] = q;
		}

		// Two passes over the neighbours: count per cell, then fill.
		int[] cellCount = new int[WORDS + 1];
		Neighbours counter = (w, from, to) -> cellCount[w + 1] += to - from;
		forEachNeighbour(query, threshold, rowMax, groupStart, grouped, counter);
		int[] start = new int[WORDS + 1];
		for (int w = 0; w < WORDS; w++) {
			start[w + 1] = start[w] + cellCount[w + 1];
		}
		int[] offsets = new int[start[WORDS]];
		int[] cursor = start.clone();
		forEachNeighbour(query, threshold, rowMax, groupStart, grouped, (w, from, to) -> {
			for (int i = from; i < to; i++) {
				offsets[cursor[w]++] = grouped[i];
			}
		});
		// One bit per word whose range is not empty: word abc is bit c of presence[ab], a table of 4 KB.
		int[] presence = new int[WORDS >>> CHAR_BITS];
		for (int w = 0; w < WORDS; w++) {
			if (start[w + 1] != start[w]) {
				presence[w >>> CHAR_BITS] |= 1 << (w & ((1 << CHAR_BITS) - 1));
			}
		}
		return new LookupTable(start, offsets, presence);
	}

	@FunctionalInterface
	private interface Neighbours {
		/** Subject word {@code w} hits the query offsets {@code grouped[from .. to)}. */
		void add(int w, int from, int to);
	}

	private static void forEachNeighbour(byte[] query, int threshold, int[] rowMax, int[] groupStart, int[] grouped,
			Neighbours out) {
		final int n = Blosum62.SIZE;
		for (int w = 0; w < WORDS; w++) {
			int from = groupStart[w], to = groupStart[w + 1];
			if (from == to) {
				continue;
			}
			int q = grouped[from];
			int a = query[q], b = query[q + 1], c = query[q + 2];
			// The query word itself is always in the table; when it scores at least the threshold the neighbour
			// enumeration below adds it.
			int self = Blosum62.score(a, a) + Blosum62.score(b, b) + Blosum62.score(c, c);
			if (self < threshold) {
				out.add(w, from, to);
			}
			// Enumerate subject words x y z, pruning with the best possible score of the remaining letters.
			int bound = rowMax[b] + rowMax[c];
			for (int x = 0; x < n; x++) {
				int sx = Blosum62.score(a, x);
				if (sx + bound < threshold) {
					continue;
				}
				for (int y = 0; y < n; y++) {
					int sxy = sx + Blosum62.score(b, y);
					if (sxy + rowMax[c] < threshold) {
						continue;
					}
					for (int z = 0; z < n; z++) {
						if (sxy + Blosum62.score(c, z) >= threshold) {
							out.add((x << (2 * CHAR_BITS)) | (y << CHAR_BITS) | z, from, to);
						}
					}
				}
			}
		}
	}

	/** Number of (word, query offset) entries. */
	public int size() {
		return offsets.length;
	}

	/** The query offsets of word {@code w}, for tests. */
	int[] offsets(int w) {
		return copyOfRange(offsets, start[w], start[w + 1]);
	}
}
