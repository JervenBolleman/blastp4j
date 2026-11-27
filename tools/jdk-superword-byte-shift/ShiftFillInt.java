/**
 * JDK 27 C2: a loop that shifts a byte and stores it into a int[] is vectorized wrongly.
 *
 *   int[] out; out[i] = in[i] << 5;
 *
 * Run: javac ShiftFillInt.java && java ShiftFillInt
 * Expect: "OK". Observed on 27+35-2325: a mismatch once the method is C2-compiled.
 * Clean under -XX:-UseSuperWord, and under -XX:TieredStopAtLevel=1.
 */
public class ShiftFillInt {

	static int[] fill(byte[] in) {
		int[] out = new int[in.length];
		for (int i = 0; i < in.length; i++) {
			out[i] = in[i] << 5;
		}
		return out;
	}

	public static void main(String[] args) {
		final int n = args.length > 0 ? Integer.parseInt(args[0]) : 142;
		final int rounds = args.length > 1 ? Integer.parseInt(args[1]) : 200_000;
		byte[] in = new byte[n];
		for (int i = 0; i < n; i++) {
			in[i] = (byte) (i % 25);
		}
		// The reference, computed before anything is hot.
		int[] want = new int[n];
		for (int i = 0; i < n; i++) {
			want[i] = in[i] << 5;
		}
		int bad = 0;
		String first = null;
		for (int r = 0; r < rounds; r++) {
			int[] got = fill(in);
			for (int i = 0; i < n; i++) {
				if (got[i] != want[i]) {
					bad++;
					if (first == null) {
						first = "round " + r + ", index " + i + ": in=" + in[i] + " want=" + want[i] + " got=" + got[i];
					}
					break;
				}
			}
		}
		System.out.println(bad == 0 ? "OK (" + rounds + " rounds, n=" + n + ")"
				: "WRONG in " + bad + " of " + rounds + " rounds, n=" + n + "; first " + first);
	}
}
