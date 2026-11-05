package swiss.sib.blastp4j.compo;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The SEG low-complexity filter of Wootton and Federhen (Comput. Chem. 17, 149 (1993)), as implemented in BLAST+'s
 * {@code blast_seg.c}, for NCBIstdaa protein sequences.
 * <p>
 * Composition-based statistics masks the subject with SEG (window 10, locut 1.8, hicut 2.1) before computing its
 * composition and redoing the alignment. The code follows BLAST's closely, including its tabulated logarithms, so that
 * the same residues are masked.
 */
public final class Seg {
	/** {@code ln(n!)} as tabulated in BLAST; beyond the table Stirling's formula is used. */
	private static final double[] LNFACT = loadLnfact();
	/** {@code ln(k / 10)}, used for the entropy of windows of exactly 10 letters. */
	private static final double[] LOG_WIN10 = { 0.0, -2.30258509, -1.60943791, -1.203982804, -0.91629073, -0.6931478,
			-0.510825623, -0.356674944, -0.22314355, -0.105360515, 0.0 };
	private static final double LN2 = 0.69314718055994530941723212145818;
	private static final double LN20 = 2.9957322735539909;
	private static final int ALPHASIZE = 20;
	private static final int FENCE_SENTRY = 201;
	/** Index in the 20-letter alphabet of each NCBIstdaa code, or -1 for letters that are not counted (X, B, Z...). */
	private static final int[] ALPHAINDEX = new int[128];

	static {
		Arrays.fill(ALPHAINDEX, -1);
		for (int c = 0, i = 0; c < 128; c++) {
			if (c == 1 || (c >= 3 && c <= 20) || c == 22) {
				ALPHAINDEX[c] = i++;
			}
		}
	}

	private final int window;
	private final double locut, hicut;
	private final int maxtrim = 50, maxbogus = 2;

	/** SEG with the given window and cutoffs; BLAST's defaults are 12, 2.2, 2.5 and CBS uses 10, 1.8, 2.1. */
	public Seg(int window, double locut, double hicut) {
		this.window = window;
		this.locut = locut;
		this.hicut = Math.max(hicut, locut);
	}

	/** The parameters composition-based statistics uses ({@code BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"}). */
	public static Seg forCompositionBasedStatistics() {
		return new Seg(10, 1.8, 2.1);
	}

	private static double[] loadLnfact() {
		List<Double> values = new ArrayList<>();
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				Seg.class.getResourceAsStream("seg_lnfact.txt"), StandardCharsets.US_ASCII))) {
			for (String line; (line = in.readLine()) != null;) {
				if (!line.isEmpty() && line.charAt(0) != '#') {
					values.add(Double.parseDouble(line));
				}
			}
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
		return values.stream().mapToDouble(Double::doubleValue).toArray();
	}

	/** A window on the sequence ({@code SSequence}). */
	private static final class Win {
		final byte[] seq;
		final Win parent;
		final int parentLength;
		int pos; // index of the first letter in seq
		int start, length, bogus;
		int[] composition, state;
		// For each count: how many classes have it, and how many classes have a larger count.
		int[] num, above;
		double entropy = -2;

		/** The whole sequence. */
		Win(byte[] seq, int pos, int length) {
			this.seq = seq;
			this.parent = null;
			this.parentLength = length;
			this.pos = pos;
			this.length = length;
		}

		/** {@code s_OpenWin}. */
		Win(Win parent, int start, int length) {
			this.seq = parent.seq;
			this.parent = parent;
			this.parentLength = parent.length;
			this.pos = parent.pos + start;
			this.start = start;
			this.length = length;
			stateOn();
		}

		void stateOn() {
			composition = new int[ALPHASIZE];
			for (int i = 0; i < length; i++) {
				int c = ALPHAINDEX[seq[pos + i]];
				if (c >= 0) {
					composition[c]++;
				} else {
					bogus++;
				}
			}
			int[] s = new int[ALPHASIZE + 1];
			int nel = 0;
			for (int c : composition) {
				if (c != 0) {
					s[nel++] = c;
				}
			}
			// descending
			Arrays.sort(s, 0, nel);
			for (int i = 0; i < nel / 2; i++) {
				int t = s[i];
				s[i] = s[nel - 1 - i];
				s[nel - 1 - i] = t;
			}
			state = s;
			num = new int[length + 2];
			above = new int[length + 2];
			for (int i = 0; i < nel; i++) {
				num[s[i]]++;
			}
			int running = 0;
			for (int v = length + 1; v >= 0; v--) {
				above[v] = running;
				running += num[v];
			}
		}

		/** One letter leaves a class that had {@code klass} letters: the last entry of that count drops by one. */
		void lower(int klass) {
			state[above[klass] + num[klass] - 1] = klass - 1;
			num[klass]--;
			num[klass - 1]++;
			above[klass - 1]--;
		}

		/** One letter joins a class that had {@code klass} letters: the first entry of that count rises by one. */
		void raise(int klass) {
			state[above[klass]] = klass + 1;
			num[klass]--;
			num[klass + 1]++;
			above[klass]++;
		}

		/** {@code s_ShiftWin1}: moves the window one letter to the right. */
		boolean shift() {
			int next = pos + length < seq.length ? seq[pos + length] & 0xff : 0;
			if (next == FENCE_SENTRY) {
				return false;
			}
			if (++start + length > parentLength) {
				--start;
				return false;
			}
			int c = ALPHAINDEX[seq[pos]];
			if (c >= 0) {
				lower(composition[c]--);
			} else {
				bogus--;
			}
			pos++;
			c = ALPHAINDEX[next];
			if (c >= 0) {
				raise(composition[c]++);
			} else {
				bogus++;
			}
			if (entropy > -2) {
				entropy = entropy(state);
			}
			return true;
		}
	}

	/** {@code s_Entropy}. */
	static double entropy(int[] sv) {
		int total = 0;
		for (int i = 0; sv[i] != 0; i++) {
			total += sv[i];
		}
		if (total == 0) {
			return 0;
		}
		double ent = 0;
		if (total == 10) {
			for (int i = 0; sv[i] != 0; i++) {
				ent += sv[i] * LOG_WIN10[sv[i]] / LN2;
			}
		} else {
			for (int i = 0; sv[i] != 0; i++) {
				ent += sv[i] * Math.log((double) sv[i] / total) / LN2;
			}
		}
		return Math.abs(ent / total);
	}

	/** {@code s_SeqEntropy}: the entropy of the window centred on each position, or -1. */
	private double[] seqEntropy(Win seq) {
		int downset = (window + 1) / 2 - 1, upset = window - downset;
		if (window > seq.length) {
			return null;
		}
		double[] h = new double[seq.length];
		Arrays.fill(h, -1);
		Win win = new Win(seq, 0, window);
		win.entropy = entropy(win.state);
		int first = downset, last = seq.length - upset;
		for (int i = first; i <= last; i++) {
			if (win.bogus > maxbogus) {
				h[i] = -1;
				win.shift();
				continue;
			}
			h[i] = win.entropy;
			win.shift();
		}
		return h;
	}

	private static int findLow(int i, int limit, double hicut, double[] h) {
		int j;
		for (j = i; j >= limit; j--) {
			if (h[j] == -1.0 || h[j] > hicut) {
				break;
			}
		}
		return j + 1;
	}

	private static int findHigh(int i, int limit, double hicut, double[] h) {
		int j;
		for (j = i; j <= limit; j++) {
			if (h[j] == -1.0 || h[j] > hicut) {
				break;
			}
		}
		return j - 1;
	}

	private static double lnfact(int n) {
		return n < LNFACT.length ? LNFACT[n] : (n + 0.5) * Math.log(n) - n + 0.9189385332;
	}

	private static double lnPerm(int[] sv, int windowLength) {
		double ans = lnfact(windowLength);
		for (int i = 0; sv[i] != 0; i++) {
			ans -= lnfact(sv[i]);
		}
		return ans;
	}

	/** {@code s_LnAss}. */
	private static double lnAss(int[] sv, int alphasize) {
		double ans = LNFACT[alphasize];
		if (sv[0] == 0) {
			return ans;
		}
		int total = alphasize, klass = 1;
		int k = 0;
		int svi = sv[0], svim1 = sv[0];
		for (int i = 0;; svim1 = svi) {
			if (++i == alphasize) {
				ans -= lnfact(klass);
				break;
			} else if ((svi = sv[++k]) == svim1) {
				klass++;
				continue;
			} else {
				total -= klass;
				ans -= lnfact(klass);
				if (svi == 0) {
					ans -= lnfact(total);
					break;
				} else {
					klass = 1;
					continue;
				}
			}
		}
		return ans;
	}

	private static double prob(int[] sv, int total) {
		double totseq = total * LN20;
		double ans1 = lnAss(sv, ALPHASIZE), ans2 = 0;
		if (ans1 > -100000.0 && sv[0] != Integer.MIN_VALUE) {
			ans2 = lnPerm(sv, total);
		}
		return ans1 + ans2 - totseq;
	}

	/** {@code s_Trim}: the sub-window of minimal probability; returns {leftend, rightend}. */
	private int[] trim(Win seq, int leftend, int rightend) {
		int lend = 0, rend = seq.length - 1;
		int minlen = 1;
		if (seq.length - maxtrim > minlen) {
			minlen = seq.length - maxtrim;
		}
		double minprob = 1.;
		for (int len = seq.length; len > minlen; len--) {
			Win win = new Win(seq, 0, len);
			boolean shift = true;
			int i = 0;
			while (shift) {
				double p = prob(win.state, len);
				if (p < minprob) {
					minprob = p;
					lend = i;
					rend = len + i - 1;
				}
				shift = win.shift();
				i++;
			}
		}
		return new int[] { leftend + lend, rightend - (seq.length - rend - 1) };
	}

	/** {@code s_SegSeq}: adds segments {begin, end} (inclusive, plus {@code offset}) to {@code segs}. */
	private void segSeq(Win seq, List<int[]> segs, int offset) {
		if (window <= 0) {
			return;
		}
		int downset = (window + 1) / 2 - 1, upset = window - downset;
		double[] h = seqEntropy(seq);
		if (h == null) {
			return;
		}
		int first = downset, last = seq.length - upset, lowlim = first;
		for (int i = first; i <= last; i++) {
			if (h[i] <= locut && h[i] != -1.0) {
				int loi = findLow(i, lowlim, hicut, h);
				int hii = findHigh(i, last, hicut, h);
				int leftend = loi - downset, rightend = hii + upset - 1;
				Win temp = new Win(seq, leftend, rightend - leftend + 1);
				int[] t = trim(temp, leftend, rightend);
				leftend = t[0];
				rightend = t[1];
				if (i + upset - 1 < leftend) {
					// a trigger window in the left trim
					int lend = loi - downset, rend = leftend - 1;
					Win left = new Win(seq, lend, rend - lend + 1);
					segSeq(left, segs, offset + lend);
				}
				segs.add(new int[] { leftend + offset, rightend + offset });
				i = Math.min(hii, rightend + downset);
				lowlim = i + 1;
			}
		}
	}

	/**
	 * The low-complexity segments of {@code seq[pos .. pos + length)}, as {begin, end} pairs (inclusive, relative to
	 * {@code pos}); they may overlap.
	 */
	public List<int[]> segments(byte[] seq, int pos, int length) {
		List<int[]> segs = new ArrayList<>();
		segSeq(new Win(seq, pos, length), segs, 0);
		return segs;
	}

	/** Replaces the low-complexity letters of {@code seq[pos .. pos + length)} by X; returns whether any were. */
	public boolean mask(byte[] seq, int pos, int length) {
		return mask(seq, pos, segments(seq, pos, length));
	}

	/** Replaces the letters of {@code segs} ({@link #segments} of {@code seq} at {@code pos}) by X; whether any are. */
	public static boolean mask(byte[] seq, int pos, List<int[]> segs) {
		for (int[] s : segs) {
			for (int i = s[0]; i <= s[1]; i++) {
				seq[pos + i] = swiss.sib.blastp4j.db.Alphabet.X;
			}
		}
		return !segs.isEmpty();
	}
}
