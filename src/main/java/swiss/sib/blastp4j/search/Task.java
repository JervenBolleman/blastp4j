package swiss.sib.blastp4j.search;

/** BLAST+'s blastp tasks that blastp4j searches with ({@code -task}). */
public enum Task {
	/** {@code blastp}: words of 3 letters, neighbourhood threshold 11. */
	BLASTP,
	/**
	 * {@code blastp-fast}: words of 5 letters of a 15-letter compressed alphabet, threshold 20
	 * ({@link swiss.sib.blastp4j.seed.CompressedTable}), and the chaining estimate before the gapped extensions
	 * ({@link Chaining}).
	 */
	BLASTP_FAST
}
