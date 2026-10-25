package swiss.sib.blastp4j.output;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import swiss.sib.blastp4j.align.EditScript;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.Alphabet;
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

	/**
	 * The columns {@code qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score
	 * qlen slen btop} for one alignment.
	 *
	 * @param query       the query (NCBIstdaa)
	 * @param residues    the array holding the subject
	 * @param subjectBase where the subject starts in {@code residues}
	 */
	public static String row(String qseqid, String sseqid, byte[] query, byte[] residues, int subjectBase,
			int subjectLength, Hsp hsp, Statistics statistics) {
		EditScript e = hsp.editScript;
		int q = hsp.queryStart, s = subjectBase + hsp.subjectStart;
		int length = 0, identities = 0, mismatches = 0, gapOpens = 0, run = 0;
		StringBuilder btop = new StringBuilder();
		for (int i = 0; i < e.size(); i++) {
			int n = e.num(i);
			length += n;
			switch (e.op(i)) {
			case EditScript.SUB -> {
				for (int k = 0; k < n; k++, q++, s++) {
					if (query[q] == residues[s]) {
						identities++;
						run++;
					} else {
						mismatches++;
						if (run > 0) {
							btop.append(run);
							run = 0;
						}
						btop.append(Alphabet.decode(query[q])).append(Alphabet.decode(residues[s]));
					}
				}
			}
			case EditScript.DEL -> { // gap in the query
				gapOpens++;
				if (run > 0) {
					btop.append(run);
					run = 0;
				}
				for (int k = 0; k < n; k++, s++) {
					btop.append('-').append(Alphabet.decode(residues[s]));
				}
			}
			default -> { // gap in the subject
				gapOpens++;
				if (run > 0) {
					btop.append(run);
					run = 0;
				}
				for (int k = 0; k < n; k++, q++) {
					btop.append(Alphabet.decode(query[q])).append('-');
				}
			}
			}
		}
		if (run > 0) {
			btop.append(run);
		}
		return String.join("\t", qseqid, sseqid,
				fixed(100.0 * identities / length, 3), Integer.toString(length),
				Integer.toString(mismatches), Integer.toString(gapOpens), Integer.toString(hsp.queryStart + 1),
				Integer.toString(hsp.queryEnd), Integer.toString(hsp.subjectStart + 1),
				Integer.toString(hsp.subjectEnd), evalue(hsp.evalue), bitScore(statistics.bitScore(hsp.score)),
				Integer.toString(hsp.score), Integer.toString(query.length), Integer.toString(subjectLength),
				btop.toString());
	}
}
