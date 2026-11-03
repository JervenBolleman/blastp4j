package swiss.sib.blastp4j.db;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * A sequence identifier ({@code Seq-id} of the NCBI data model), as BLAST databases store them and as BLAST+ prints
 * them.
 * <p>
 * Identifiers come from the {@code .phr} file of a BLAST database, or from a FASTA header parsed as
 * {@code makeblastdb -parse_seqids} does. BLAST prints them in three ways, reproduced here (checked against blastp
 * for each kind):
 * <ul>
 * <li>{@code sseqid} ({@link #seqid}): all identifiers in FASTA form, joined with {@code |}, for example
 * {@code sp|P69905|HBA_HUMAN} or {@code gi|123|ref|NP_000509.1|}; local identifiers without {@code lcl|};</li>
 * <li>{@code sacc} and {@code saccver} ({@link #accession}): the accession of the best identifier, without or with
 * its version: {@code P69905}, {@code NP_000509}, {@code 1ABC_A} for a PDB chain, {@code db:tag} for a general
 * identifier.</li>
 * </ul>
 *
 * @param kind      the kind of identifier
 * @param accession accession, local string, PDB molecule or general database
 * @param name      locus name, PDB chain or general tag; may be null
 * @param number    version of an accession, or a number (gi, numeric local id); -1 if none
 * @param release   release ({@code reviewed} or {@code unreviewed} for UniProt); may be null
 */
public record SeqId(Kind kind, String accession, String name, long number, String release) {
	/** The choices of {@code Seq-id}, in ASN.1 tag order. */
	public enum Kind {
		LOCAL, GIBBSQ, GIBBMT, GIIM, GENBANK, EMBL, PIR, SWISSPROT, PATENT, OTHER, GENERAL, GI, DDBJ, PRF, PDB, TPG, TPE,
		TPD, GPIPE, NAMED_ANNOT_TRACK
	}

	private static boolean isTextseq(Kind k) {
		return switch (k) {
		case GENBANK, EMBL, PIR, SWISSPROT, OTHER, DDBJ, PRF, TPG, TPE, TPD, GPIPE, NAMED_ANNOT_TRACK -> true;
		default -> false;
		};
	}

	/** The FASTA tag of a kind of text identifier. */
	private String tag() {
		return switch (kind) {
		case GENBANK -> "gb";
		case EMBL -> "emb";
		case PIR -> "pir";
		case SWISSPROT -> "unreviewed".equals(release) ? "tr" : "sp";
		case OTHER -> "ref";
		case DDBJ -> "dbj";
		case PRF -> "prf";
		case TPG -> "tpg";
		case TPE -> "tpe";
		case TPD -> "tpd";
		case GPIPE -> "gpp";
		case NAMED_ANNOT_TRACK -> "nat";
		default -> throw new IllegalStateException(kind.toString());
		};
	}

	private String accessionVersion() {
		return number > 0 ? accession + "." + number : accession;
	}

	/** The identifier in FASTA form, as part of {@code sseqid}. */
	public String fasta() {
		if (isTextseq(kind)) {
			return tag() + "|" + (accession == null ? "" : accessionVersion()) + "|" + (name == null ? "" : name);
		}
		return switch (kind) {
		case LOCAL -> accession != null ? accession : Long.toString(number);
		case GI -> "gi|" + number;
		case GIBBSQ -> "bbs|" + number;
		case GIBBMT -> "bbm|" + number;
		case PDB -> "pdb|" + accession + "|" + (name == null ? "" : name);
		case GENERAL -> "gnl|" + accession + "|" + name;
		default -> accession != null ? accession : Long.toString(number);
		};
	}

	/** The accession, with or without version. */
	public String accession(boolean withVersion) {
		if (isTextseq(kind)) {
			return withVersion ? accessionVersion() : accession;
		}
		return switch (kind) {
		case PDB -> name == null || name.isEmpty() ? accession : accession + "_" + name;
		case GENERAL -> accession + ":" + name;
		case LOCAL -> accession != null ? accession : Long.toString(number);
		default -> Long.toString(number);
		};
	}

	/** How good an identifier is for an accession: text identifiers and PDB first, then general, local, gi. */
	private int rank() {
		if (isTextseq(kind) || kind == Kind.PDB) {
			return 0;
		}
		return switch (kind) {
		case GENERAL -> 1;
		case LOCAL -> 2;
		default -> 3;
		};
	}

	/** {@code sseqid}: all identifiers in FASTA form. */
	public static String seqid(List<SeqId> ids) {
		StringBuilder sb = new StringBuilder();
		for (SeqId id : ids) {
			if (!sb.isEmpty()) {
				sb.append('|');
			}
			sb.append(id.fasta());
		}
		return sb.toString();
	}

	/** {@code sacc} or {@code saccver}: the accession of the best identifier. */
	public static String accession(List<SeqId> ids, boolean withVersion) {
		SeqId best = ids.get(0);
		for (SeqId id : ids) {
			if (id.rank() < best.rank()) {
				best = id;
			}
		}
		return best.accession(withVersion);
	}

	/**
	 * Parses the first word of a FASTA header as {@code makeblastdb -parse_seqids} does: a sequence of NCBI FASTA
	 * identifiers ({@code sp|P69905|HBA_HUMAN}, {@code gi|123|ref|NP_000509.1|}, {@code lcl|x}, ...); a word that is
	 * not in that syntax is a local identifier.
	 */
	public static List<SeqId> parseFasta(String word) {
		List<SeqId> ids = new ArrayList<>();
		String[] f = word.split("\\|", -1);
		int i = 0;
		while (i < f.length) {
			String tag = f[i].toLowerCase(Locale.ROOT);
			SeqId id = null;
			int fields = 0;
			Kind text = switch (tag) {
			case "gb" -> Kind.GENBANK;
			case "emb" -> Kind.EMBL;
			case "pir" -> Kind.PIR;
			case "sp", "tr" -> Kind.SWISSPROT;
			case "ref" -> Kind.OTHER;
			case "dbj" -> Kind.DDBJ;
			case "prf" -> Kind.PRF;
			case "tpg" -> Kind.TPG;
			case "tpe" -> Kind.TPE;
			case "tpd" -> Kind.TPD;
			case "gpp" -> Kind.GPIPE;
			case "nat" -> Kind.NAMED_ANNOT_TRACK;
			default -> null;
			};
			if (text != null && i + 1 < f.length) {
				String acc = f[i + 1];
				long version = -1;
				int dot = acc.lastIndexOf('.');
				if (dot > 0 && dot < acc.length() - 1 && acc.substring(dot + 1).chars().allMatch(Character::isDigit)) {
					version = Long.parseLong(acc.substring(dot + 1));
					acc = acc.substring(0, dot);
				}
				String name = i + 2 < f.length && !f[i + 2].isEmpty() ? f[i + 2] : null;
				String release = text == Kind.SWISSPROT ? (tag.equals("sp") ? "reviewed" : "unreviewed") : null;
				id = new SeqId(text, acc, name, version, release);
				fields = i + 2 < f.length ? 3 : 2;
			} else if ((tag.equals("gi") || tag.equals("bbs") || tag.equals("bbm")) && i + 1 < f.length
					&& isNumber(f[i + 1])) {
				Kind k = tag.equals("gi") ? Kind.GI : tag.equals("bbs") ? Kind.GIBBSQ : Kind.GIBBMT;
				id = new SeqId(k, null, null, Long.parseLong(f[i + 1]), null);
				fields = 2;
			} else if (tag.equals("lcl") && i + 1 < f.length) {
				id = new SeqId(Kind.LOCAL, f[i + 1], null, -1, null);
				fields = 2;
			} else if (tag.equals("gnl") && i + 2 < f.length) {
				id = new SeqId(Kind.GENERAL, f[i + 1], f[i + 2], -1, null);
				fields = 3;
			} else if (tag.equals("pdb") && i + 1 < f.length) {
				String chain = i + 2 < f.length ? f[i + 2] : null;
				id = new SeqId(Kind.PDB, f[i + 1], chain, -1, null);
				fields = i + 2 < f.length ? 3 : 2;
			}
			if (id == null) {
				// Not NCBI FASTA syntax: the whole word is a local identifier.
				return List.of(new SeqId(Kind.LOCAL, word, null, -1, null));
			}
			ids.add(id);
			i += fields;
		}
		return ids.isEmpty() ? List.of(new SeqId(Kind.LOCAL, word, null, -1, null)) : ids;
	}

	private static boolean isNumber(String s) {
		return !s.isEmpty() && s.chars().allMatch(Character::isDigit);
	}

	/** An identifier from its ASN.1 {@code Seq-id} choice. */
	static SeqId fromAsn(Asn1.Node choice) {
		Kind kind = Kind.values()[choice.tag()];
		Asn1.Node v = choice.inner();
		if (isTextseq(kind)) {
			// Textseq-id ::= SEQUENCE { name [0], accession [1], release [2], version [3] }
			return new SeqId(kind, str(v.child(1)), str(v.child(0)), v.child(3) == null ? -1 : v.child(3).inner().integer(),
					str(v.child(2)));
		}
		return switch (kind) {
		case LOCAL -> objectId(kind, null, v);
		case GENERAL -> {
			// Dbtag ::= SEQUENCE { db [0] VisibleString, tag [1] Object-id }
			Asn1.Node tag = v.child(1).inner();
			String db = str(v.child(0));
			yield tag.tag() == 1 ? new SeqId(kind, db, tag.inner().string(), -1, null)
					: new SeqId(kind, db, Long.toString(tag.inner().integer()), -1, null);
		}
		case PDB -> {
			// PDB-seq-id ::= SEQUENCE { mol [0], chain [1] INTEGER, rel [2], chain-id [3] VisibleString }
			String chain = v.child(3) != null ? str(v.child(3))
					: v.child(1) != null ? String.valueOf((char) v.child(1).inner().integer()) : null;
			yield new SeqId(kind, str(v.child(0)), chain, -1, null);
		}
		default -> new SeqId(kind, null, null, v.constructed() ? -1 : v.integer(), null);
		};
	}

	/** {@code Object-id ::= CHOICE { id [0] INTEGER, str [1] VisibleString }} */
	private static SeqId objectId(Kind kind, String db, Asn1.Node choice) {
		return choice.tag() == 1 ? new SeqId(kind, choice.inner().string(), null, -1, null)
				: new SeqId(kind, null, null, choice.inner().integer(), null);
	}

	private static String str(Asn1.Node tagged) {
		return tagged == null ? null : tagged.inner().string();
	}

	@Override
	public String toString() {
		return fasta();
	}
}
