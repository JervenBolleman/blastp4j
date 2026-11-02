package swiss.sib.blastp4j.stats;

import swiss.sib.blastp4j.db.Alphabet;

/**
 * The ungapped Karlin-Altschul parameters lambda, K and H of a query, as blastp computes them
 * ({@code Blast_ScoreBlkKbpUngappedCalc} in BLAST+'s {@code blast_stat.c}). They depend on the query's residue
 * composition: the score distribution is that of BLOSUM62 between the query's letters (X not counted) and the
 * Robinson and Robinson background frequencies.
 * <p>
 * blastp prints them as the first "Lambda K H" line of its text report. It uses them for the ungapped stage: the X-drop
 * of the ungapped extensions and the gap trigger, the score an ungapped alignment needs to be extended with gaps.
 *
 * @param lambda the scale of the scores
 * @param k      the K parameter
 * @param h      the relative entropy, in nats
 */
public record UngappedKarlin(double lambda, double k, double h) {
	/** Background frequencies of Robinson and Robinson (1991), in NCBIstdaa order of ACDEFGHIKLMNPQRSTVWY. */
	private static final String ROBINSON_LETTERS = "ACDEFGHIKLMNPQRSTVWY";
	private static final double[] ROBINSON = { 78.05, 19.25, 53.64, 62.95, 38.56, 73.77, 21.99, 51.42, 57.44, 90.19,
			22.43, 44.87, 52.03, 42.64, 51.29, 71.20, 58.41, 64.41, 13.30, 32.16 };

	private static final double LAMBDA_ACCURACY = 1e-5;
	private static final int LAMBDA_ITERATIONS = 17;
	private static final double LAMBDA_START = 0.5;
	private static final double K_SUM_LIMIT = 0.0001;
	private static final int K_ITERATIONS = 100;

	/**
	 * The parameters of BLOSUM62 with the Robinson and Robinson composition on both sides ({@code kbp_ideal} of
	 * BLAST).
	 */
	public static UngappedKarlin ideal() {
		double[] background = robinson();
		return of(new ScoreDistribution(background, background));
	}

	/** The Robinson and Robinson background frequencies per NCBIstdaa code, summing to 1 (0 for the other codes). */
	public static double[] robinson() {
		double[] background = new double[Blosum62.SIZE];
		for (int i = 0; i < ROBINSON.length; i++) {
			background[Alphabet.encode(ROBINSON_LETTERS.charAt(i))] = ROBINSON[i];
		}
		normalize(background);
		return background;
	}

	public double logK() {
		return Math.log(k);
	}

	/** The parameters of a query given in NCBIstdaa codes. */
	public static UngappedKarlin of(byte[] query) {
		return of(query, 0, query.length);
	}

	/** The parameters of {@code length} NCBIstdaa codes of {@code query} from {@code offset}. */
	public static UngappedKarlin of(byte[] query, int offset, int length) {
		int n = Blosum62.SIZE;
		double[] queryFreq = new double[n];
		for (int i = offset; i < offset + length; i++) {
			queryFreq[query[i]]++;
		}
		queryFreq[Alphabet.X] = 0; // ambiguous residues are not counted
		normalize(queryFreq);
		return of(new ScoreDistribution(queryFreq, robinson()));
	}

	private static void normalize(double[] p) {
		double sum = 0;
		for (double v : p) {
			sum += v;
		}
		if (sum > 0) {
			for (int i = 0; i < p.length; i++) {
				p[i] /= sum;
			}
		}
	}

	static UngappedKarlin of(ScoreDistribution d) {
		if (d.average >= 0) {
			throw new IllegalArgumentException("the expected score must be negative");
		}
		double lambda = lambda(d);
		double h = entropy(d, lambda);
		double k = k(d, lambda, h);
		return new UngappedKarlin(lambda, k, h);
	}

	/**
	 * Probabilities of each pair score, the {@code BlastScoreFreqCalc} of BLAST: {@code prob[score - min]} for scores
	 * from {@code min} to {@code max}, where only scores with non-zero probability count for the bounds.
	 */
	static final class ScoreDistribution {
		final int lo, min, max;
		final double[] prob; // indexed by score - lo
		final double average;

		ScoreDistribution(double[] p1, double[] p2) {
			int lo = Integer.MAX_VALUE, hi = Integer.MIN_VALUE;
			for (int a = 0; a < Blosum62.SIZE; a++) {
				for (int b = 0; b < Blosum62.SIZE; b++) {
					lo = Math.min(lo, Blosum62.score(a, b));
					hi = Math.max(hi, Blosum62.score(a, b));
				}
			}
			this.lo = lo;
			prob = new double[hi - lo + 1];
			for (int a = 0; a < Blosum62.SIZE; a++) {
				for (int b = 0; b < Blosum62.SIZE; b++) {
					prob[Blosum62.score(a, b) - lo] += p1[a] * p2[b];
				}
			}
			double sum = 0;
			int min = Integer.MIN_VALUE, max = Integer.MIN_VALUE;
			for (int s = lo; s <= hi; s++) {
				if (prob[s - lo] > 0) {
					sum += prob[s - lo];
					max = s;
					if (min == Integer.MIN_VALUE) {
						min = s;
					}
				}
			}
			this.min = min;
			this.max = max;
			double avg = 0;
			if (sum > 0.0001 || sum < -0.0001) {
				for (int s = min; s <= max; s++) {
					prob[s - lo] /= sum;
					avg += s * prob[s - lo];
				}
			}
			this.average = avg;
		}

		double p(int score) {
			return prob[score - lo];
		}
	}

	private static int gcd(int a, int b) {
		a = Math.abs(a);
		b = Math.abs(b);
		while (b != 0) {
			int t = a % b;
			a = b;
			b = t;
		}
		return a;
	}

	/** {@code Blast_KarlinLambdaNR}: the positive root of sum_s p(s) e^(lambda s) = 1, by safeguarded Newton. */
	static double lambda(ScoreDistribution d) {
		return lambda(d.prob, d.lo, d.min, d.max, d.average, LAMBDA_START);
	}

	/**
	 * {@code Blast_KarlinLambdaNR} for score probabilities {@code prob[s - offset]} between {@code low} and {@code high}
	 * with average score {@code average}, starting from {@code lambda0}; -1 when lambda does not exist.
	 */
	public static double lambda(double[] prob, int offset, int low, int high, double average, double lambda0) {
		if (average >= 0 || low >= 0 || high <= 0) {
			return -1;
		}
		int g = -low;
		for (int i = 1; i <= high - low && g > 1; i++) {
			if (prob[i + low - offset] != 0) {
				g = gcd(g, i);
			}
		}
		final int itmax = 20, maxNewton = 20 + LAMBDA_ITERATIONS;
		double a = 0, b = 1, f = 4;
		double x0 = Math.exp(-lambda0);
		double x = 0 < x0 && x0 < 1 ? x0 : 0.5;
		boolean isNewton = false;
		for (int k = 0; k < itmax; k++) {
			double fold = f;
			boolean wasNewton = isNewton;
			isNewton = false;
			// Horner's rule for the polynomial and its derivative
			double gd = 0;
			f = prob[low - offset];
			for (int i = low + g; i < 0; i += g) {
				gd = x * gd + f;
				f = f * x + prob[i - offset];
			}
			gd = x * gd + f;
			f = f * x + prob[-offset] - 1;
			for (int i = g; i <= high; i += g) {
				gd = x * gd + f;
				f = f * x + prob[i - offset];
			}
			if (f > 0) {
				a = x;
			} else if (f < 0) {
				b = x;
			} else {
				break;
			}
			if (b - a < 2 * a * (1 - b) * LAMBDA_ACCURACY) {
				x = (a + b) / 2;
				break;
			}
			if (k >= maxNewton || (wasNewton && Math.abs(f) > .9 * Math.abs(fold)) || gd >= 0) {
				x = (a + b) / 2;
			} else {
				double p = -f / gd;
				double y = x + p;
				if (y <= a || y >= b) {
					x = (a + b) / 2;
				} else {
					isNewton = true;
					x = y;
					if (Math.abs(p) < LAMBDA_ACCURACY * x * (1 - x)) {
						break;
					}
				}
			}
		}
		return -Math.log(x) / g;
	}

	/** {@code BlastKarlinLtoH}: the relative entropy. */
	static double entropy(ScoreDistribution d, double lambda) {
		int low = d.min, high = d.max;
		double etonlam = Math.exp(-lambda);
		double sum = low * d.p(low);
		for (int s = low + 1; s <= high; s++) {
			sum = s * d.p(s) + etonlam * sum;
		}
		double scale = Math.pow(etonlam, high);
		return scale > 0 ? lambda * sum / scale : lambda * Math.exp(lambda * high + Math.log(sum));
	}

	/** {@code BlastKarlinLHtoK}: K from lambda and H, by the series of Karlin and Altschul. */
	static double k(ScoreDistribution d, double lambda, double h) {
		int low = d.min, high = d.max;
		int range = high - low;
		int divisor = -low;
		for (int i = 1; i <= range && divisor > 1; i++) {
			if (d.p(low + i) != 0) {
				divisor = gcd(divisor, i);
			}
		}
		high /= divisor;
		low /= divisor;
		lambda *= divisor;
		range = high - low;
		double firstTermClosedForm = h / lambda;
		double expMinusLambda = Math.exp(-lambda);
		if (low == -1 && high == 1) {
			double pl = d.p(low * divisor), ph = d.p(high * divisor);
			return (pl - ph) * (pl - ph) / pl;
		}
		if (low == -1 || high == 1) {
			if (high != 1) {
				double avg = d.average / divisor;
				firstTermClosedForm = avg * avg / firstTermClosedForm;
			}
			return firstTermClosedForm * (1 - expMinusLambda);
		}
		// probArrayStartLow[i] = probability of score (low + i) * divisor... as in BLAST, indexed from the lowest
		// observed score without dividing (BLAST indexes sprob[low] after low was divided; for divisor 1 it is equal).
		double[] probs = new double[range + 1];
		for (int i = 0; i <= range; i++) {
			int s = d.min + i;
			probs[i] = s <= d.max ? d.p(s) : 0;
		}
		double[] asp = new double[K_ITERATIONS * range + 1];
		double outerSum = 0, innerSum = 1;
		int lowAlignmentScore = 0, highAlignmentScore = 0;
		asp[0] = 1;
		for (int iter = 0; iter < K_ITERATIONS && innerSum > K_SUM_LIMIT;) {
			int first = range, last = range;
			lowAlignmentScore += low;
			highAlignmentScore += high;
			int p;
			for (p = highAlignmentScore - lowAlignmentScore; p >= 0; p--) {
				int p1 = p - first, p1e = p - last, p2 = first;
				innerSum = 0;
				for (; p1 >= p1e; p1--, p2++) {
					innerSum += asp[p1] * probs[p2];
				}
				if (first != 0) {
					--first;
				}
				if (p <= range) {
					--last;
				}
				asp[p] = innerSum;
			}
			p++; // the loop left p at -1; BLAST's pointer is incremented before use
			innerSum = asp[p];
			int i;
			for (i = lowAlignmentScore + 1; i < 0; i++) {
				innerSum = asp[++p] + innerSum * expMinusLambda;
			}
			innerSum *= expMinusLambda;
			for (; i <= highAlignmentScore; ++i) {
				innerSum += asp[++p];
			}
			innerSum /= ++iter;
			outerSum += innerSum;
		}
		return -Math.exp(-2.0 * outerSum) / (firstTermClosedForm * Math.expm1(-lambda));
	}
}
