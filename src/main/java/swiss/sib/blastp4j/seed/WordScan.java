package swiss.sib.blastp4j.seed;

import java.util.List;

/**
 * A word finder's scans, for a {@link WordFinder} or a copy of it from {@link WordFinder#of}: the caller calls
 * {@link #findBlocks} when {@link #vectorScan} says so, otherwise {@link #find}.
 */
public interface WordScan {
	/** Whether {@link #findBlocks} is the faster scan for this query. */
	boolean vectorScan();

	/** {@link WordFinder#find}. */
	void find(byte[] residues, int start, int length, List<UngappedHsp> out);

	/** {@link WordFinder#findBlocks}. */
	void findBlocks(byte[] residues, int start, int length, List<UngappedHsp> out);
}
