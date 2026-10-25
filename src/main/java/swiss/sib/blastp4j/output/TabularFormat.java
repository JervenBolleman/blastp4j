package swiss.sib.blastp4j.output;

import java.util.Locale;

/** Number formatting of BLAST+ tabular output ({@code -outfmt 6}). */
public final class TabularFormat {
	private TabularFormat() {
	}

	/** E-value: {@code 0.0}, {@code 1.23e-45}, {@code 0.012}, {@code 0.12}, {@code 1.2} or a whole number. */
	public static String evalue(double evalue) {
		if (evalue < 1.0e-180) {
			return "0.0";
		}
		if (evalue < 0.0009) {
			return String.format(Locale.ROOT, "%.2e", evalue);
		}
		if (evalue < 0.1) {
			return String.format(Locale.ROOT, "%.3f", evalue);
		}
		if (evalue < 1.0) {
			return String.format(Locale.ROOT, "%.2f", evalue);
		}
		if (evalue < 10.0) {
			return String.format(Locale.ROOT, "%.1f", evalue);
		}
		return String.format(Locale.ROOT, "%.0f", evalue);
	}

	/** Bit score: one decimal up to 99.9, then the integer part, and exponent notation above 9999. */
	public static String bitScore(double bitScore) {
		if (bitScore > 9999) {
			return String.format(Locale.ROOT, "%.3e", bitScore);
		}
		if (bitScore > 99.9) {
			return Long.toString((long) bitScore);
		}
		return String.format(Locale.ROOT, "%.1f", bitScore);
	}
}
