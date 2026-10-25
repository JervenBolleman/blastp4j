package swiss.sib.blastp4j.stats;

/**
 * The complementary error function, which Java does not provide. Relative error about 1e-15 over the whole range.
 * <p>
 * For |x| below 1.2 it is {@code 1 - erf(x)} with erf from its Taylor series (no cancellation, erfc is above 0.08
 * there). Otherwise {@code erfc(x) = Q(1/2, x^2)}, the regularised upper incomplete gamma function, evaluated with its
 * continued fraction (modified Lentz). Negative arguments use {@code erfc(-x) = 2 - erfc(x)}.
 */
final class Erfc {
	private static final double SQRT_PI = Math.sqrt(Math.PI);
	private static final double EPS = 1e-16;
	private static final double TINY = 1e-300;

	private Erfc() {
	}

	static double erfc(double x) {
		if (x < 0) {
			return 2.0 - erfc(-x);
		}
		if (x < 1.2) {
			return 1.0 - erfSeries(x);
		}
		return upperGammaHalf(x * x);
	}

	/** erf(x) = 2/sqrt(pi) * sum (-1)^n x^(2n+1) / (n! (2n+1)). */
	private static double erfSeries(double x) {
		double x2 = x * x;
		double term = x; // (-1)^n x^(2n+1) / n!
		double sum = x;
		for (int n = 1; n < 100; n++) {
			term *= -x2 / n;
			double add = term / (2 * n + 1);
			sum += add;
			if (Math.abs(add) < EPS * Math.abs(sum)) {
				break;
			}
		}
		return 2.0 / SQRT_PI * sum;
	}

	/** Q(1/2, z) by the continued fraction of the incomplete gamma function. */
	private static double upperGammaHalf(double z) {
		final double a = 0.5;
		double b = z + 1.0 - a;
		double c = 1.0 / TINY;
		double d = 1.0 / b;
		double h = d;
		for (int i = 1; i < 10_000; i++) {
			double an = -i * (i - a);
			b += 2.0;
			d = an * d + b;
			if (Math.abs(d) < TINY) {
				d = TINY;
			}
			c = b + an / c;
			if (Math.abs(c) < TINY) {
				c = TINY;
			}
			d = 1.0 / d;
			double del = d * c;
			h *= del;
			if (Math.abs(del - 1.0) < EPS) {
				break;
			}
		}
		// exp(-z) z^a / Gamma(a) with Gamma(1/2) = sqrt(pi)
		return Math.exp(-z + a * Math.log(z)) / SQRT_PI * h;
	}
}
