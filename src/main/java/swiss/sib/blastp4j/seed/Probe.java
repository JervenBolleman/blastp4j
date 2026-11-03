package swiss.sib.blastp4j.seed;

import java.util.List;

/** A query's part of a batch's scan: {@link WordFinder#probe} and {@link WordFinder#endSubject}, or a copy's. */
interface Probe {
	void probe(byte[] residues, int start, int length, int q, int s, List<UngappedHsp> out);

	void endSubject(int length);
}
