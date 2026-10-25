package swiss.sib.blastp4j.db;

/**
 * The identifiers and titles of the sequences of a database, decoded when asked for: from FASTA headers or from the
 * {@code .phr} file of a BLAST database. Implementations are thread safe.
 */
public interface Headers {
	/** Number of sequences. */
	int size();

	/** The first word of the header, as written in the FASTA file. */
	String name(int i);

	/** {@code sseqid}: the identifiers as BLAST prints them. */
	String seqid(int i);

	/** {@code sacc} (without version) or {@code saccver}. */
	String accession(int i, boolean withVersion);

	/** {@code stitle}: the title (description). */
	String title(int i);
}
