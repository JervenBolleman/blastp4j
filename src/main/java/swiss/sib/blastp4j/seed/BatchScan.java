package swiss.sib.blastp4j.seed;


/** The scan of a batch of queries, for a {@link BatchWordFinder} or a copy of it from {@link BatchWordFinder#of}. */
public interface BatchScan {
	/** {@link BatchWordFinder#find}. */
	void find(byte[] residues, int start, int length, SeedList[] out);
}
