package swiss.sib.blastp4j.util;

/**
 * How the floating-point code computes: as BLAST+ does, operation for operation, so that the results are blastp's to
 * the last bit; or faster, with fused multiply-adds and sums in another order, which can change the last bits and so,
 * rarely, a result. blastp4j uses {@link #BLAST}, almostblastp {@link #FAST}.
 */
public enum Arithmetic {
	/** BLAST+'s operations in BLAST+'s order. */
	BLAST,
	/** Fused multiply-adds and reordered sums where they are faster. */
	FAST
}
