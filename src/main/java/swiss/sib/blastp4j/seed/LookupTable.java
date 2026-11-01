package swiss.sib.blastp4j.seed;

import static java.util.Arrays.copyOfRange;
import swiss.sib.blastp4j.stats.Blosum62;

/**
 * blastp's lookup table: for every word of {@code wordSize} letters (three for blastp), the query positions whose word
 * scores at least the threshold against it (the neighbourhood words), plus the query's own words whatever their score
 * ({@code blast_aalookup.c}).
 * <p>
 * A word is indexed by its NCBIstdaa codes, five bits each, first letter highest, as in BLAST. The table is stored
 * compressed: the query offsets of word {@code w} are {@code offsets[start[w] .. start[w + 1])}.
 * <p>
 * Instances are immutable.
 */
public record LookupTable(int wordSize, int[] start, int[] offsets) {
	/** Letters per word. */
	public static final int WORD_SIZE = 3;
	/** Bits per letter. */
	static final int CHAR_BITS = 5;
	/** Number of words of {@link #WORD_SIZE} letters (some unused, the alphabet has 28 letters). */
	static final int WORDS = 1 << (WORD_SIZE * CHAR_BITS);

	/** Index of the word of {@code wordSize} letters at {@code pos} of {@code seq}. */
	static int word(byte[] seq, int pos, int wordSize) {
		int word = 0;
		for (int i = 0; i < wordSize; i++) {
			word = (word << CHAR_BITS) | seq[pos + i];
		}
		return word;
	}

	/** The table of a query (NCBIstdaa codes) for neighbourhood score {@code threshold} (11 for blastp), with three letters per word. */
	public static LookupTable of(byte[] query, int threshold) {
		return of(query, threshold, WORD_SIZE);
	}

	/** The table of a query (NCBIstdaa codes) for neighbourhood score {@code threshold} and {@code wordSize} letters per word. */
	public static LookupTable of(byte[] query, int threshold, int wordSize) {
		final int n = Blosum62.SIZE;
		final int cells = 1 << (wordSize * CHAR_BITS);
		int[] rowMax = new int[n];
		for (int a = 0; a < n; a++) {
			rowMax[a] = Integer.MIN_VALUE;
			for (int b = 0; b < n; b++) {
				rowMax[a] = Math.max(rowMax[a], Blosum62.score(a, b));
			}
		}
		// Group the query offsets by their word, then find the neighbours once per distinct word.
		int words = Math.max(0, query.length - wordSize + 1);
		int[] count = new int[cells + 1];
		for (int q = 0; q < words; q++) {
			count[word(query, q, wordSize) + 1]++;
		}
		int[] groupStart = new int[cells + 1];
		for (int w = 0; w < cells; w++) {
			groupStart[w + 1] = groupStart[w] + count[w + 1];
		}
		int[] grouped = new int[words];
		int[] fill = groupStart.clone();
		for (int q = 0; q < words; q++) {
			grouped[fill[word(query, q, wordSize)]++] = q;
		}

		// Two passes over the neighbours: count per cell, then fill.
		int[] cellCount = new int[cells + 1];
		Neighbours counter = (w, from, to) -> cellCount[w + 1] += to - from;
		forEachNeighbour(query, threshold, wordSize, rowMax, groupStart, grouped, counter);
		int[] start = new int[cells + 1];
		for (int w = 0; w < cells; w++) {
			start[w + 1] = start[w] + cellCount[w + 1];
		}
		int[] offsets = new int[start[cells]];
		int[] cursor = start.clone();
		forEachNeighbour(query, threshold, wordSize, rowMax, groupStart, grouped, (w, from, to) -> {
			for (int i = from; i < to; i++) {
				offsets[cursor[w]++] = grouped[i];
			}
		});
		return new LookupTable(wordSize, start, offsets);
	}

	@FunctionalInterface
	private interface Neighbours {
		/** Subject word {@code w} hits the query offsets {@code grouped[from .. to)}. */
		void add(int w, int from, int to);
	}

	private static void forEachNeighbour(byte[] query, int threshold, int wordSize, int[] rowMax, int[] groupStart,
			int[] grouped, Neighbours out) {
		final int cells = groupStart.length - 1;
		// bound[i]: the best possible score of the query letters from i on, to prune the enumeration.
		final int[] bound = new int[wordSize + 1];
		for (int w = 0; w < cells; w++) {
			int from = groupStart[w], to = groupStart[w + 1];
			if (from == to) {
				continue;
			}
			int q = grouped[from];
			// The query word itself is always in the table; when it scores at least the threshold the neighbour
			// enumeration below adds it.
			int self = 0;
			for (int i = 0; i < wordSize; i++) {
				self += Blosum62.score(query[q + i], query[q + i]);
			}
			if (self < threshold) {
				out.add(w, from, to);
			}
			for (int i = wordSize - 1; i >= 0; i--) {
				bound[i] = bound[i + 1] + rowMax[query[q + i]];
			}
			enumerate(query, q, 0, 0, 0, wordSize, threshold, bound, from, to, out);
		}
	}

	/** Enumerates the subject words that score at least the threshold against the query word at {@code q}. */
	private static void enumerate(byte[] query, int q, int position, int score, int code, int wordSize, int threshold,
			int[] bound, int from, int to, Neighbours out) {
		if (position == wordSize) {
			out.add(code, from, to);
			return;
		}
		final int shift = (wordSize - 1 - position) * CHAR_BITS;
		for (int x = 0; x < Blosum62.SIZE; x++) {
			int sx = score + Blosum62.score(query[q + position], x);
			if (sx + bound[position + 1] >= threshold) {
				enumerate(query, q, position + 1, sx, code | (x << shift), wordSize, threshold, bound, from, to, out);
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
