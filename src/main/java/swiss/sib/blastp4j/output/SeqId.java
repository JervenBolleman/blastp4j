package swiss.sib.blastp4j.output;

import java.util.Set;

/**
 * Accessions of sequence identifiers in NCBI FASTA syntax, as BLAST+ prints them for a database built with
 * {@code makeblastdb -parse_seqids} ({@code sacc}, {@code saccver}).
 * <p>
 * Queries are read without {@code -parse_deflines}, so BLAST treats their identifiers as local ones and prints them
 * whole for {@code qacc} and {@code qaccver}.
 */
public final class SeqId {
	/** Identifier types whose second field is an accession, possibly with a version. */
	private static final Set<String> ACCESSION_TYPES = Set.of("sp", "tr", "ref", "gb", "emb", "dbj", "pir", "prf",
			"tpg", "tpe", "tpd", "gpp", "nat");

	private SeqId() {
	}

	/**
	 * The accession of an identifier: {@code sp|P69905|HBA_HUMAN} gives {@code P69905}, {@code ref|NP_000509.1|}
	 * gives {@code NP_000509.1} (or {@code NP_000509} without version), {@code pdb|1ABC|A} gives {@code 1ABC_A} and
	 * {@code lcl|x} gives {@code x}. Other identifiers are returned unchanged.
	 */
	public static String accession(String id, boolean withVersion) {
		String[] f = id.split("\\|", -1);
		if (f.length < 2) {
			return id;
		}
		String type = f[0].toLowerCase(java.util.Locale.ROOT);
		String acc;
		if (ACCESSION_TYPES.contains(type)) {
			acc = f[1];
		} else if (type.equals("pdb") && f.length >= 3 && !f[2].isEmpty()) {
			return f[1] + "_" + f[2];
		} else if (type.equals("lcl") || type.equals("pdb")) {
			return f[1];
		} else {
			return id;
		}
		if (!withVersion) {
			int dot = acc.lastIndexOf('.');
			if (dot > 0) {
				acc = acc.substring(0, dot);
			}
		}
		return acc;
	}
}
