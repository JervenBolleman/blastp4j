package swiss.sib.blastp4j.align;

/**
 * Many independent {@link GappedAligner#alignEx} extensions, submitted and then advanced a row at a time, each with
 * exactly alignEx's result: {@link TracebackLanes} (in vector lanes) or {@link TracebackPairs} (two at a time, in
 * scalar code). Not thread safe: use one per thread.
 */
public interface TracebackKernel {
	/** Whether this takes extensions with these gap costs and A of up to {@code maxA} letters. */
	boolean fits(int gapOpen, int gapExtend, int maxA);

	/** The longest A an extension may have. */
	int maxA();

	/** Extensions submitted and not yet started. */
	int queued();

	/**
	 * Submits {@code alignEx(a, aBase, m, b, bBase, n, reverse, xDrop)}, with {@code 0 < m <= maxA} and {@code n > 0},
	 * and the score matrix {@code matrix[a << 5 | b]} (32 x 32) of the given {@code version}: the matrix is copied
	 * only where its version is not the one last copied, so the caller may refill one array with another matrix under
	 * a new version, once no extension it submitted with the old one is queued. {@link #step} reports the result
	 * under {@code id}.
	 */
	void submit(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, int id,
			int[] matrix, long version);

	/**
	 * Advances every extension by one row, reporting those that end to {@code done} (which may submit more). Returns
	 * false, doing nothing, when there is no extension, started or queued.
	 */
	boolean step(TracebackLanes.Done done);
}
