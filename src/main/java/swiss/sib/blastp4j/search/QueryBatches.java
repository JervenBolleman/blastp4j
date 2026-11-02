package swiss.sib.blastp4j.search;

/**
 * BLAST's query batches, which decide where the composition-based stage may skip subjects.
 * <p>
 * blastp searches its queries in batches: it adds queries to a batch until the batch holds at least 10,000 residues
 * ({@code GetQueryBatchSize}, {@code CBlastInput::GetNextSeqBatch}; the {@code BATCH_SIZE} environment variable
 * overrides the size). With composition-based statistics, a subject is skipped only when the heaps of all queries in
 * the batch are filled to their cutoff ({@code BlastCompo_EarlyTermination}). With one thread the queries are done in
 * order, so this can only happen for the last query of a batch, when all earlier ones ended filled.
 * <p>
 * Usage: call {@link #next} for each query in order, search with its answer, then call {@link #done} with whether
 * the search's heap ended filled ({@link BlastpSearch.Result#filledToCutoff()}).
 */
public final class QueryBatches {
	/** blastp's batch size in residues. */
	public static final int BLASTP_BATCH_SIZE = 10000;

	private final int[] lengths;
	private final int batchSize;
	private int next;
	private long batchResidues;
	private boolean earlierFilled = true;

	/**
	 * @param lengths   the query lengths, in order
	 * @param batchSize the batch size in residues
	 */
	public QueryBatches(int[] lengths, int batchSize) {
		this.lengths = lengths.clone();
		this.batchSize = batchSize;
	}

	/** Batches of {@link #BLASTP_BATCH_SIZE} residues, or of {@code BATCH_SIZE} from the environment as in BLAST. */
	public static QueryBatches blastp(int[] lengths) {
		String size = System.getenv("BATCH_SIZE");
		return new QueryBatches(lengths, size != null ? Integer.parseInt(size.trim()) : BLASTP_BATCH_SIZE);
	}

	/** The batch size in residues: a batch takes queries until it holds at least this many. */
	public int batchSize() {
		return batchSize;
	}

	/** Whether the next query may skip subjects: it is the last of its batch, and all earlier ones ended filled. */
	public boolean next() {
		batchResidues += lengths[next];
		boolean last = batchResidues >= batchSize || next + 1 == lengths.length;
		return last && earlierFilled;
	}

	/** Records the end of the current query's search: whether its heap ended filled to the cutoff. */
	public void done(boolean filledToCutoff) {
		earlierFilled &= filledToCutoff;
		if (batchResidues >= batchSize) {
			batchResidues = 0;
			earlierFilled = true;
		}
		next++;
	}
}
