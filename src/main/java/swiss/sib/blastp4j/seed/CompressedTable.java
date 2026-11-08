package swiss.sib.blastp4j.seed;

import java.util.Arrays;

/**
 * The lookup table of {@code blastp -task blastp-fast} for a batch of queries: BLAST's compressed-alphabet table
 * ({@code BlastCompressedAaLookupTableNew} in {@code blast_aalookup.c}) over the queries concatenated, as BLAST builds
 * it for a batch. A word is 5 letters of the {@link CompressedAlphabet}, so the table has 15^5 cells; word
 * {@code c0 c1 c2 c3 c4} is cell {@code c0 + 15 c1 + 225 c2 + 3375 c3 + 50625 c4}.
 * <p>
 * Each query word without a letter outside the alphabet adds its offset to the cells of its neighbours: the compressed
 * words that score at least the threshold against it ({@link CompressedAlphabet#SCORES}, scaled by 100, as is the
 * threshold), and to its own compressed word if that scores less.
 * <p>
 * The hits of cell {@code w} are entries {@code [cells[2 w], cells[2 w + 2])}, in the order BLAST's scan reads them: a cell
 * keeps its first 5 offsets in the backbone and the rest in cells of 4, the newest first, so after the first 3 come
 * the others in chunks of 4, the last chunk first ({@code s_CompressedLookupAddWordHit},
 * {@code s_BlastCompressedAaScanSubject}). That order is the order of a subject position's hits, which decides ties
 * later. An entry is the offset in the concatenation, {@code base[i]} plus the query offset; {@code query[k]} is the
 * query of entry {@code k}. {@code cells[2 w + 1]} is the cell's first entry again (0 where it has none): next to
 * its start, so that a position whose cell has one hit (most have one or two) reads one cache line of the table, not
 * two, as BLAST's backbone cell holds its first hits. {@code size} is the length of the diagonal array.
 * <p>
 * {@code narrow} holds the same entries as chars, half the memory, when the concatenation is shorter than 65,536
 * letters (BLAST's batches are about 10,000), and is null otherwise; the scan then reads it, as
 * {@link BatchWordFinder.Table#narrow} for blastp.
 * <p>
 * {@code present} has a bit per cell, set where the cell has hits ({@code present[w >>> 6] >>> w & 1}), so that the
 * scan reads 95 KB for a position without hits rather than the 6 MB of {@code cells}, as BLAST's scan tests its
 * {@code pv} array: half the cells are empty, and two thirds of the positions (bean queries) fall on one.
 */
public record CompressedTable(int[] cells, int[] entry, char[] narrow, int[] query, int[] base, int size,
		long[] present) {
	/** Word length. */
	public static final int WORD = 5;
	/** Cells: 15^5. */
	static final int CELLS = 759375;
	/** {@code BLAST_WORD_THRESHOLD_BLASTP_FAST}. */
	public static final int THRESHOLD = 20;
	/** Hits a backbone cell holds before it overflows, and that stay first in the scan. */
	private static final int BACKBONE = 5, FIRST = 3, CHUNK = 4;

	public static CompressedTable of(byte[][] queries, int threshold, int window) {
		final int n = queries.length;
		final int[] base = new int[n];
		long length = 0;
		for (int i = 0; i < n; i++) {
			if (i > 0 && queries[i - 1].length > 0) {
				length++;
			}
			base[i] = (int) length;
			length += queries[i].length;
		}
		if (length + window > 1 << 30) {
			throw new IllegalArgumentException("a batch of " + n + " queries is too long: " + length);
		}
		final Neighbours words = new Neighbours((int) (CompressedAlphabet.SCALE * threshold));
		for (int i = 0; i < n; i++) {
			words.query = i;
			for (int q = 0; q + WORD <= queries[i].length; q++) {
				words.add(queries[i], q, base[i] + q);
			}
		}
		// The cells in insertion order, then each in the scan's order.
		final int[] start = new int[CELLS + 1];
		for (int h = 0; h < words.count; h++) {
			start[words.cell[h] + 1]++;
		}
		for (int c = 0; c < CELLS; c++) {
			start[c + 1] += start[c];
		}
		final int[] fill = Arrays.copyOf(start, CELLS), entry = new int[words.count], query = new int[words.count];
		for (int h = 0; h < words.count; h++) {
			final int k = fill[words.cell[h]]++;
			entry[k] = words.offset[h];
			query[k] = words.queryOf[h];
		}
		final int[] e = new int[CHUNK * 64], qy = new int[CHUNK * 64];
		final long[] present = new long[(CELLS + 63) >>> 6];
		for (int c = 0; c < CELLS; c++) {
			final int from = start[c], count = start[c + 1] - from;
			if (count > 0) {
				present[c >>> 6] |= 1L << c;
			}
			if (count > BACKBONE) {
				reorder(entry, query, from, count, e, qy);
			}
		}
		final int[] cells = new int[2 * CELLS + 2];
		for (int c = 0; c < CELLS; c++) {
			cells[2 * c] = start[c];
			cells[2 * c + 1] = start[c + 1] > start[c] ? entry[start[c]] : 0;
		}
		cells[2 * CELLS] = start[CELLS];
		char[] narrow = null;
		if (length <= Character.MAX_VALUE) {
			narrow = new char[entry.length];
			for (int j = 0; j < entry.length; j++) {
				narrow[j] = (char) entry[j];
			}
		}
		return new CompressedTable(cells, entry, narrow, query, base, WordFinder.diagonals((int) length, window),
				present);
	}

	/** A cell's entries after the first 3, in chunks of 4 in insertion order: the chunks reversed. */
	private static void reorder(int[] entry, int[] query, int from, int count, int[] e, int[] qy) {
		final int rest = count - FIRST, first = from + FIRST;
		if (e.length < rest) {
			e = new int[rest];
			qy = new int[rest];
		}
		System.arraycopy(entry, first, e, 0, rest);
		System.arraycopy(query, first, qy, 0, rest);
		int to = first;
		for (int chunk = (rest - 1) / CHUNK * CHUNK; chunk >= 0; chunk -= CHUNK) {
			final int len = Math.min(CHUNK, rest - chunk);
			System.arraycopy(e, chunk, entry, to, len);
			System.arraycopy(qy, chunk, query, to, len);
			to += len;
		}
	}

	/** The neighbouring words of query words ({@code s_CompressedAddWordHits} and its core), as (cell, offset) pairs. */
	private static final class Neighbours {
		final int threshold;
		/** Per amino acid code, the best score against any compressed letter. */
		final int[] rowMax = new int[CompressedAlphabet.COMPRESS.length];
		/** Per code, the scores against the compressed letters from best to worst, and those letters. */
		final int[][] sorted = new int[rowMax.length][], letters = new int[rowMax.length][];
		int[] cell = new int[1024], offset = new int[1024], queryOf = new int[1024];
		int count, query;
		private final int[] word = new int[WORD], subject = new int[WORD];
		private int at;

		Neighbours(int threshold) {
			this.threshold = threshold;
			final int size = CompressedAlphabet.SIZE;
			for (int a = 0; a < rowMax.length; a++) {
				int top = CompressedAlphabet.SCORES[a * size];
				for (int j = 1; j < size; j++) {
					top = Math.max(top, CompressedAlphabet.SCORES[a * size + j]);
				}
				rowMax[a] = top;
				final int max = top;
				// s_loadSortedMatrix: by the difference from the row's best, stable.
				final Integer[] order = new Integer[size];
				for (int j = 0; j < size; j++) {
					order[j] = j;
				}
				final int row = a;
				Arrays.sort(order, (x, y) -> Integer.compare(max - CompressedAlphabet.SCORES[row * size + x],
						max - CompressedAlphabet.SCORES[row * size + y]));
				sorted[a] = new int[size];
				letters[a] = new int[size];
				for (int j = 0; j < size; j++) {
					letters[a][j] = order[j];
					sorted[a][j] = CompressedAlphabet.SCORES[a * size + order[j]];
				}
			}
		}

		/** The query word at {@code q}, at offset {@code at} in the concatenation. */
		void add(byte[] query, int q, int at) {
			int self = 0;
			for (int i = 0; i < WORD; i++) {
				final int c = CompressedAlphabet.COMPRESS[query[q + i]];
				if (c >= CompressedAlphabet.SIZE) {
					return; // a letter outside the alphabet: the word is not indexed
				}
				word[i] = query[q + i];
				subject[i] = c;
				self += CompressedAlphabet.SCORES[query[q + i] * CompressedAlphabet.SIZE + c];
			}
			this.at = at;
			if (self < threshold) {
				save(index(subject));
			}
			int bound = 0;
			for (int i = 0; i < WORD; i++) {
				bound += rowMax[word[i]];
			}
			neighbours(bound, 0);
		}

		/** {@code s_CompressedAddWordHitsCore}: {@code score} is the exact score so far plus the best of the rest. */
		private void neighbours(int score, int pos) {
			final int a = word[pos];
			score -= rowMax[a];
			final int[] row = sorted[a], letter = letters[a];
			for (int i = 0; i < CompressedAlphabet.SIZE && score + row[i] >= threshold; i++) {
				subject[pos] = letter[i];
				if (pos == WORD - 1) {
					save(index(subject));
				} else {
					neighbours(score + row[i], pos + 1);
				}
			}
		}

		private static int index(int[] w) {
			return w[0] + 15 * (w[1] + 15 * (w[2] + 15 * (w[3] + 15 * w[4])));
		}

		private void save(int c) {
			if (count == cell.length) {
				cell = Arrays.copyOf(cell, 2 * count);
				offset = Arrays.copyOf(offset, 2 * count);
				queryOf = Arrays.copyOf(queryOf, 2 * count);
			}
			cell[count] = c;
			offset[count] = at;
			queryOf[count++] = query;
		}
	}
}
