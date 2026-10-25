package swiss.sib.blastp4j.output;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

import swiss.sib.blastp4j.align.EditScript;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.stats.Blosum62;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * BLAST+ tabular output ({@code -outfmt 6}): its number formatting and the columns of one alignment.
 * <p>
 * Numbers are rounded as C's {@code printf} does: the exact binary value of the double, ties to even. Java's
 * {@code String.format} rounds ties up, which differs for values such as 25/64 = 39.0625 ("39.062" in BLAST).
 */
public final class TabularFormat {
	private TabularFormat() {
	}

	/** {@code printf("%.<decimals>f", x)}. */
	static String fixed(double x, int decimals) {
		return new BigDecimal(x).setScale(decimals, RoundingMode.HALF_EVEN).toPlainString();
	}

	/** {@code printf("%.<decimals>e", x)} for x &gt; 0. */
	static String exponent(double x, int decimals) {
		BigDecimal d = new BigDecimal(x).round(new MathContext(decimals + 1, RoundingMode.HALF_EVEN));
		int exp = d.precision() - d.scale() - 1;
		BigDecimal mantissa = d.movePointLeft(exp).setScale(decimals, RoundingMode.HALF_EVEN);
		return mantissa.toPlainString() + "e" + (exp < 0 ? "-" : "+") + (Math.abs(exp) < 10 ? "0" : "")
				+ Math.abs(exp);
	}

	/** E-value: {@code 0.0}, {@code 1.23e-45}, {@code 0.012}, {@code 0.12}, {@code 1.2} or a whole number. */
	public static String evalue(double evalue) {
		if (evalue < 1.0e-180) {
			return "0.0";
		}
		if (evalue < 0.0009) {
			return exponent(evalue, 2);
		}
		if (evalue < 0.1) {
			return fixed(evalue, 3);
		}
		if (evalue < 1.0) {
			return fixed(evalue, 2);
		}
		if (evalue < 10.0) {
			return fixed(evalue, 1);
		}
		return fixed(evalue, 0);
	}

	/** Bit score: one decimal up to 99.9, then the integer part, and exponent notation above 9999. */
	public static String bitScore(double bitScore) {
		if (bitScore > 9999) {
			return exponent(bitScore, 3);
		}
		if (bitScore > 99.9) {
			return Long.toString((long) bitScore);
		}
		return fixed(bitScore, 1);
	}

	/** The columns of {@code -outfmt 6} that are supported, with BLAST's names. */
	public enum Column {
		QSEQID, SSEQID, QACC, QACCVER, SACC, SACCVER, QLEN, SLEN, QSTART, QEND, SSTART, SEND, QSEQ, SSEQ, EVALUE, BITSCORE, SCORE, LENGTH, PIDENT,
		NIDENT, MISMATCH, POSITIVE, GAPOPEN, GAPS, PPOS, QCOVHSP, BTOP, STITLE;

		/** The column for a BLAST field name, or null. */
		public static Column of(String name) {
			for (Column c : values()) {
				if (c.name().equalsIgnoreCase(name)) {
					return c;
				}
			}
			return null;
		}
	}

	/** {@code -outfmt 6} without field names ("std"). */
	public static final List<Column> STANDARD = List.of(Column.QACCVER, Column.SACCVER, Column.PIDENT, Column.LENGTH,
			Column.MISMATCH, Column.GAPOPEN, Column.QSTART, Column.QEND, Column.SSTART, Column.SEND, Column.EVALUE,
			Column.BITSCORE);

	/**
	 * Parses the field list of {@code -outfmt "6 ..."}; "std" stands for the standard columns.
	 *
	 * @throws IllegalArgumentException for an unknown or unsupported field
	 */
	public static List<Column> parseColumns(String fields) {
		List<Column> columns = new ArrayList<>();
		for (String f : fields.trim().split("\\s+")) {
			if (f.isEmpty()) {
				continue;
			}
			if (f.equals("std")) {
				columns.addAll(STANDARD);
				continue;
			}
			Column c = Column.of(f);
			if (c == null) {
				throw new IllegalArgumentException("unsupported output field '" + f + "'; supported: std, "
						+ Arrays.stream(Column.values()).map(x -> x.name().toLowerCase(Locale.ROOT))
								.collect(Collectors.joining(" ")));
			}
			columns.add(c);
		}
		return columns.isEmpty() ? STANDARD : columns;
	}

	/**
	 * Query coverage in percent, as BLAST rounds it: to the nearest integer below 99%, truncated from there, so that
	 * 100 means complete coverage.
	 */
	static int coverage(int covered, int length) {
		double c = 100.0 * covered / length;
		return (int) (c < 99 ? c + 0.5 : c);
	}

	/** Counts and strings of one alignment, computed from its edit script. */
	private static final class Alignment {
		int length, identities, mismatches, positives, gapOpens, gaps;
		final StringBuilder btop = new StringBuilder(), qseq = new StringBuilder(), sseq = new StringBuilder();

		Alignment(byte[] query, byte[] residues, int subjectBase, Hsp hsp, boolean sequences) {
			EditScript e = hsp.editScript;
			int q = hsp.queryStart, s = subjectBase + hsp.subjectStart, run = 0;
			for (int i = 0; i < e.size(); i++) {
				int n = e.num(i);
				length += n;
				if (e.op(i) == EditScript.SUB) {
					for (int k = 0; k < n; k++, q++, s++) {
						char a = Alphabet.decode(query[q]), b = Alphabet.decode(residues[s]);
						if (sequences) {
							qseq.append(a);
							sseq.append(b);
						}
						if (Blosum62.score(query[q], residues[s]) > 0) {
							positives++;
						}
						if (query[q] == residues[s]) {
							identities++;
							run++;
						} else {
							mismatches++;
							if (run > 0) {
								btop.append(run);
								run = 0;
							}
							btop.append(a).append(b);
						}
					}
				} else {
					gapOpens++;
					gaps += n;
					if (run > 0) {
						btop.append(run);
						run = 0;
					}
					boolean inQuery = e.op(i) == EditScript.DEL;
					for (int k = 0; k < n; k++) {
						char a = inQuery ? '-' : Alphabet.decode(query[q++]);
						char b = inQuery ? Alphabet.decode(residues[s++]) : '-';
						btop.append(a).append(b);
						if (sequences) {
							qseq.append(a);
							sseq.append(b);
						}
					}
				}
			}
			if (run > 0) {
				btop.append(run);
			}
		}
	}

	/**
	 * One line of {@code -outfmt 6} for one alignment of the query with subject {@code subject} of {@code db}.
	 *
	 * @param qseqid the query's id (the first word of its header)
	 * @param query  the query (NCBIstdaa)
	 */
	public static String row(List<Column> columns, String qseqid, byte[] query, SequenceDatabase db, int subject,
			Hsp hsp, Statistics statistics) {
		final byte[] residues = db.residues();
		final int subjectBase = db.start(subject), subjectLength = db.length(subject);
		boolean sequences = columns.contains(Column.QSEQ) || columns.contains(Column.SSEQ);
		Alignment a = new Alignment(query, residues, subjectBase, hsp, sequences);
		StringBuilder line = new StringBuilder();
		for (Column c : columns) {
			if (!line.isEmpty()) {
				line.append('\t');
			}
			line.append(switch (c) {
			case QSEQID, QACC, QACCVER -> qseqid;
			case SSEQID -> db.id(subject);
			case SACC -> db.accession(subject, false);
			case SACCVER -> db.accession(subject, true);
			case STITLE -> db.title(subject);
			case QLEN -> Integer.toString(query.length);
			case SLEN -> Integer.toString(subjectLength);
			case QSTART -> Integer.toString(hsp.queryStart + 1);
			case QEND -> Integer.toString(hsp.queryEnd);
			case SSTART -> Integer.toString(hsp.subjectStart + 1);
			case SEND -> Integer.toString(hsp.subjectEnd);
			case QSEQ -> a.qseq;
			case SSEQ -> a.sseq;
			case EVALUE -> evalue(hsp.evalue);
			case BITSCORE -> bitScore(statistics.bitScore(hsp.score));
			case SCORE -> Integer.toString(hsp.score);
			case LENGTH -> Integer.toString(a.length);
			case PIDENT -> fixed(100.0 * a.identities / a.length, 3);
			case NIDENT -> Integer.toString(a.identities);
			case MISMATCH -> Integer.toString(a.mismatches);
			case POSITIVE -> Integer.toString(a.positives);
			case GAPOPEN -> Integer.toString(a.gapOpens);
			case GAPS -> Integer.toString(a.gaps);
			case PPOS -> fixed(100.0 * a.positives / a.length, 2);
			case QCOVHSP -> Integer.toString(coverage(hsp.queryEnd - hsp.queryStart, query.length));
			case BTOP -> a.btop;
			});
		}
		return line.toString();
	}

	/** The columns of the reference files: {@code qseqid sseqid pident length mismatch gapopen qstart qend sstart send
	 * evalue bitscore score qlen slen btop}. */
	public static final List<Column> REFERENCE = List.of(Column.QSEQID, Column.SSEQID, Column.PIDENT, Column.LENGTH,
			Column.MISMATCH, Column.GAPOPEN, Column.QSTART, Column.QEND, Column.SSTART, Column.SEND, Column.EVALUE,
			Column.BITSCORE, Column.SCORE, Column.QLEN, Column.SLEN, Column.BTOP);

	/** {@link #row(List, String, byte[], SequenceDatabase, int, Hsp, Statistics)} with the reference columns. */
	public static String row(String qseqid, byte[] query, SequenceDatabase db, int subject, Hsp hsp,
			Statistics statistics) {
		return row(REFERENCE, qseqid, query, db, subject, hsp, statistics);
	}

}
