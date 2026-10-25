package swiss.sib.blastp4j;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.zip.GZIPInputStream;

/** The BLAST+ reference results in src/test/resources/reference (see reference/README.md). */
public final class ReferenceData {
	/** Swiss-Prot 2026_03 as built by makeblastdb: total residues and sequences. */
	public static final long SPROT_RESIDUES = 209_017_843L;
	public static final int SPROT_SEQUENCES = 575_748;

	/** Column indexes of the tabular reference files. */
	public static final int QSEQID = 0, SSEQID = 1, EVALUE = 10, BITSCORE = 11, SCORE = 12, QLEN = 13, SLEN = 14;

	private ReferenceData() {
	}

	/** Rows of a reference file, split into its 16 columns. */
	public static List<String[]> load(String name) throws IOException {
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(ReferenceData.class.getResourceAsStream("/reference/" + name)),
				StandardCharsets.US_ASCII))) {
			return in.lines().map(l -> l.split("\t", -1)).toList();
		}
	}
}
