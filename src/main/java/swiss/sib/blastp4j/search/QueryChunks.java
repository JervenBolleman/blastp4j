package swiss.sib.blastp4j.search;

/**
 * BLAST's splitting of a query batch into chunks ({@code CQuerySplitter} and {@code SplitQuery_CalculateNumChunks} in
 * {@code split_query_cxx.cpp} and {@code split_query_aux_priv.cpp}).
 * <p>
 * BLAST runs the preliminary stage once per chunk, each chunk being a piece of the batch's concatenated queries, and
 * then maps the hits back to the whole query and merges them ({@code BlastHSPStreamMerge}). This is not only a
 * division of labour: a chunk is a query of its own to the preliminary stage, so its seeding and its cutoff scores are
 * those of a query of the chunk's length, which are lower than the whole query's. A long query therefore keeps
 * alignments that the whole query's cutoff would drop, which is why blastp4j, searching the query in one piece,
 * reported a different set of HSPs for titin (34,350 residues) until it split as well.
 * <p>
 * The split happens only when the batch holds at least {@code 2 * (CHUNK_SIZE - OVERLAP)} residues, 19,800 by default,
 * so it is reached by a single query longer than that (one is alone in its batch) or by a batch whose last query takes
 * it past the size. Chunks overlap by {@link #OVERLAP} residues so that an alignment across a boundary is found whole
 * in one of them.
 */
public final class QueryChunks {
	/** BLAST's chunk size in residues for a protein search ({@code SplitQuery_GetChunkSize}). */
	public static final int CHUNK_SIZE = 10000;
	/** BLAST's overlap between chunks ({@code SplitQuery_GetOverlapChunkSize}; {@code OVERLAP_CHUNK_SIZE} overrides). */
	public static final int OVERLAP = 100;
	/**
	 * Largest difference between the end diagonal of one chunk's HSP and the start diagonal of the next chunk's for
	 * the two to be merged ({@code OVERLAP_DIAG_CLOSE}).
	 */
	public static final int MERGE_DIAG_CLOSE = 10;

	private final int[] starts, ends;
	private final int overlap;

	private QueryChunks(int[] starts, int[] ends, int overlap) {
		this.starts = starts;
		this.ends = ends;
		this.overlap = overlap;
	}

	/** The chunks of a batch of {@code total} residues; one chunk, the whole batch, when BLAST does not split it. */
	public static QueryChunks of(int total) {
		return of(total, CHUNK_SIZE, overlapSize());
	}

	/** {@code OVERLAP_CHUNK_SIZE} from the environment as in BLAST, or {@link #OVERLAP}. */
	static int overlapSize() {
		String size = System.getenv("OVERLAP_CHUNK_SIZE");
		return size != null && !size.isBlank() ? Integer.parseInt(size.trim()) : OVERLAP;
	}

	/**
	 * The chunks of {@code total} residues for a chunk size and an overlap.
	 *
	 * @param total     the residues of the batch's concatenated queries
	 * @param chunkSize the chunk size before BLAST evens it out ({@link #CHUNK_SIZE})
	 * @param overlap   the overlap between chunks ({@link #OVERLAP})
	 */
	public static QueryChunks of(int total, int chunkSize, int overlap) {
		// SplitQuery_CalculateNumChunks: the number of chunks is the count of non-overlapping pieces, and with fewer
		// than two of them there is no split. Integer division, so 19,799 residues is one chunk and 19,800 is two.
		int chunks = chunkSize > overlap ? total / (chunkSize - overlap) : 0;
		if (chunks <= 1) {
			return new QueryChunks(new int[] { 0 }, new int[] { total }, overlap);
		}
		// The size is then evened out over the chunks, and rounded up where that keeps their number.
		int size = (total + (chunks - 1) * overlap) / chunks;
		if (chunks < size - overlap) {
			size++;
		}
		// x_ComputeChunkRanges: each chunk starts a size minus an overlap after the one before, and the last one, or
		// any one that would run past the end, ends at the end.
		int[] starts = new int[chunks], ends = new int[chunks];
		int count = 0;
		for (int start = 0, chunk = 0; chunk < chunks; chunk++) {
			int end = start + size;
			if (end >= total || chunk + 1 == chunks) {
				end = total;
			}
			starts[count] = start;
			ends[count] = end;
			count++;
			start += size - overlap;
			if (start > total || end == total) {
				break;
			}
		}
		return new QueryChunks(java.util.Arrays.copyOf(starts, count), java.util.Arrays.copyOf(ends, count), overlap);
	}

	/** The number of chunks; 1 when the batch is not split. */
	public int count() {
		return starts.length;
	}

	/** Whether BLAST splits this batch, that is whether there is more than one chunk. */
	public boolean split() {
		return starts.length > 1;
	}

	/** The first residue of a chunk, as an offset in the concatenated queries. */
	public int start(int chunk) {
		return starts[chunk];
	}

	/** One past the last residue of a chunk. */
	public int end(int chunk) {
		return ends[chunk];
	}

	/** The overlap these chunks were made with. */
	public int overlap() {
		return overlap;
	}

	@Override
	public String toString() {
		StringBuilder s = new StringBuilder("QueryChunks[");
		for (int i = 0; i < starts.length; i++) {
			s.append(i > 0 ? ", " : "").append(starts[i]).append("..").append(ends[i]);
		}
		return s.append(']').toString();
	}
}
