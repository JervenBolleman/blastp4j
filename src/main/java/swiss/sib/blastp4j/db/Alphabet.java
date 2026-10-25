package swiss.sib.blastp4j.db;

/**
 * BLAST's protein alphabet NCBIstdaa: the codes used in BLAST databases ({@code .psq}) and in its scoring matrices.
 * Code 0 is the gap character, which BLAST also uses as the sentinel between sequences.
 */
public final class Alphabet {
	/** The NCBIstdaa letters, indexed by code. */
	public static final String LETTERS = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
	/** Number of codes. */
	public static final int SIZE = LETTERS.length();
	/** The sentinel / gap code. */
	public static final byte GAP = 0;
	/** The code of X, used for letters outside the alphabet. */
	public static final byte X = (byte) LETTERS.indexOf('X');

	private static final byte[] ENCODE = new byte[128];

	static {
		java.util.Arrays.fill(ENCODE, X);
		for (int code = 0; code < SIZE; code++) {
			char c = LETTERS.charAt(code);
			ENCODE[c] = (byte) code;
			ENCODE[Character.toLowerCase(c)] = (byte) code;
		}
	}

	private Alphabet() {
	}

	/** The code of a letter; upper and lower case are the same, anything else is X. */
	public static byte encode(int letter) {
		return letter >= 0 && letter < 128 ? ENCODE[letter] : X;
	}

	/** The letter of a code. */
	public static char decode(int code) {
		return LETTERS.charAt(code);
	}

	/** Encodes a whole sequence. */
	public static byte[] encode(CharSequence sequence) {
		byte[] out = new byte[sequence.length()];
		for (int i = 0; i < out.length; i++) {
			out[i] = encode(sequence.charAt(i));
		}
		return out;
	}
}
