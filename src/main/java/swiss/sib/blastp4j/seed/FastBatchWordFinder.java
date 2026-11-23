package swiss.sib.blastp4j.seed;

import java.util.Arrays;
import java.util.List;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShape;
import jdk.incubator.vector.VectorSpecies;

/**
 * The word finder of {@code blastp -task blastp-fast} for a batch of queries: {@link BatchWordFinder} with BLAST's
 * compressed-alphabet scan ({@code s_BlastCompressedAaScanSubject}) and its {@link CompressedTable}. A subject position
 * whose 5 letters are all in the {@link CompressedAlphabet} is looked up; the others (with {@code X}, {@code U},
 * {@code *}, {@code O} or {@code -}) are skipped, as BLAST skips them. The two-hit logic is {@link BatchWordFinder}'s
 * with words of 5 letters, and the hits that are not just stored go to the query's {@link FastWordFinder}.
 * <p>
 * Not thread safe: use one per thread; the {@link CompressedTable} can be shared.
 */
public final class FastBatchWordFinder implements BatchScan {
	private static final int WORD = CompressedTable.WORD;
	/** The weight of a word's newest letter: 15^4. */
	private static final int NEWEST = 50625;
	/** How many positions with hits ahead the scan touches the cells of. */
	private static final int AHEAD = 16;

	/** The machine's preferred vector of ints, one subject position per lane, as {@link WordFinder}'s. */
	private static final VectorSpecies<Integer> INTS = IntVector.SPECIES_PREFERRED;
	/** Positions per vector. */
	private static final int LANES = INTS.length();
	/** The letters of {@link #LANES} positions before they are widened: the smallest byte vector with that many. */
	private static final VectorSpecies<Byte> BYTES = VectorSpecies.of(byte.class,
			VectorShape.forBitSize(Math.max(64, LANES * Byte.SIZE)));
	/**
	 * Whether {@link #positionsWithHits} scans with vectors. Not in a GraalVM native image, where a Vector API gather
	 * compiles to scalar loads (as in {@link WordFinder}).
	 */
	private static final boolean VECTOR = System.getProperty("org.graalvm.nativeimage.imagecode") == null;
	/** Positions per block, one bit each in a long. A multiple of {@link #LANES}, so only the last block is partial. */
	private static final int BLOCK = 64;
	/**
	 * The compressed letter of a code outside the alphabet: a bit that no compressed letter has, so that the five
	 * letters of a word are tested together with one {@code and}, and larger than any of them, so that a word with one
	 * of them is outside the table.
	 */
	private static final int OUTSIDE = 16;
	/**
	 * {@link CompressedAlphabet#COMPRESS} with {@link #OUTSIDE} for the codes without a compressed letter, padded to a
	 * power of two so that C2 drops the bounds check of the lookup.
	 */
	private static final byte[] LETTERS = new byte[32];
	static {
		Arrays.fill(LETTERS, (byte) OUTSIDE);
		for (int code = 0; code < CompressedAlphabet.COMPRESS.length; code++) {
			final byte letter = CompressedAlphabet.COMPRESS[code];
			if (letter < CompressedAlphabet.SIZE) {
				LETTERS[code] = letter;
			}
		}
	}

	private final CompressedTable table;
	private final FastWordFinder[] finders;
	private final int[] diags;
	private final int window, diagMask;
	private int diagOffset;
	/** The positions with hits of the subject being scanned, and their words. */
	private int[] words = new int[0], at = new int[0];
	/** The subject's letters compressed, with room for the loads of the last vector past its end. */
	private byte[] letters = new byte[0];
	/** One block's words. */
	private final int[] blockWords = new int[BLOCK + LANES];
	/** The sum of the cells touched ahead, kept so that the loads stay. */
	private int sink;

	public FastBatchWordFinder(byte[][] queries, SeedParameters[] parameters, CompressedTable table) {
		this.table = table;
		window = parameters[0].window();
		diags = new int[table.size()];
		diagMask = table.size() - 1;
		finders = new FastWordFinder[queries.length];
		for (int i = 0; i < queries.length; i++) {
			if (parameters[i].window() != window) {
				throw new IllegalArgumentException("the queries of a batch have different windows");
			}
			finders[i] = new FastWordFinder(queries[i], parameters[i], diags, table.base()[i]);
		}
		diagOffset = window;
	}

	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp>[] out) {
		final int[] cells = table.cells(), entry = table.entry(), query = table.query(), base = table.base();
		final int[] present = table.present();
		final int[] d = diags;
		final int diagMask = this.diagMask & (d.length - 1);
		final int offset = diagOffset, windowL = window, lastL = window - WORD - 1;
		final FastWordFinder[] f = finders;
		// First the positions whose cell has hits, with their words; then their hits, each position touching the cells
		// of the one AHEAD further on, so that the misses on cells (6 MB) overlap rather than wait one after the other.
		if (words.length < length) {
			words = new int[length];
			at = new int[length];
		}
		final int[] w = words, pos = at;
		final int n = VECTOR ? positionsWithHits(residues, start, length, present)
				: scalarPositionsWithHits(residues, start, length, present);
		int touched = 0;
		for (int j = 0; j < n; j++) {
			if (j + AHEAD < n) {
				touched += cells[2 * w[j + AHEAD]];
			}
			final int cell = 2 * w[j], s = pos[j];
			// The cell's first entry from its own line of cells, the others from entry. A counted loop, which C2
			// keeps without a safepoint poll per hit: each hit loads the next one's entry (entry has one more).
			int e = cells[cell + 1];
			for (int k = cells[cell], end = cells[cell + 2]; k < end; k++) {
				// As BatchWordFinder.find.
				final int diag = (e - s) & diagMask;
				final int state = d[diag];
				final int diff = s - (state - offset), x = diff - WORD;
				if ((~(x | (lastL - x)) | state << 1) >= 0) {
					final int older = (diff - windowL) >> 31;
					d[diag] = (state & older) | ((s + offset) & ~older);
				} else {
					final int i = query[k];
					f[i].hit(residues, start, length, e - base[i], s, out[i]);
				}
				e = entry[k + 1];
			}
		}
		sink += touched;
		for (FastWordFinder finder : f) {
			finder.endSubject(length);
		}
		diagOffset = diagOffset >= Integer.MAX_VALUE / 4 ? window : diagOffset + length + window;
	}

	/**
	 * The starts of the words of {@code residues[start .. start + length)} whose cell has hits, in {@link #at}, and
	 * their words in {@link #words}; how many there are. BLAST's scan ({@code s_BlastCompressedAaScanSubject}) rolls
	 * the word one position at a time and branches on the presence bit, which is unpredictable: two thirds of the
	 * positions have no hits, in no pattern. Here a block of 64 positions gives its words and presence bits without
	 * branching, {@link #LANES} at a time, and only the positions of the block's mask are written out.
	 */
	private int positionsWithHits(byte[] residues, int start, int length, int[] present) {
		final byte[] c = compressed(residues, start, length);
		final int[] w = words, pos = at, bw = blockWords;
		final int positions = length - WORD + 1;
		int n = 0;
		for (int block = 0; block < positions; block += BLOCK) {
			long mask = blockMask(c, block, Math.min(BLOCK, positions - block), present, bw);
			for (; mask != 0; mask &= mask - 1) {
				final int b = Long.numberOfTrailingZeros(mask);
				w[n] = bw[b];
				pos[n++] = block + b;
			}
		}
		return n;
	}

	/**
	 * The subject compressed, so that a word is five loads and no table lookup. The letters past its end are outside
	 * the alphabet, which makes the words that would run past it outside the table.
	 */
	private byte[] compressed(byte[] residues, int start, int length) {
		if (letters.length < length + BLOCK) {
			letters = new byte[length + BLOCK];
		}
		final byte[] c = letters;
		for (int p = 0; p < length; p++) {
			c[p] = LETTERS[residues[start + p] & 31];
		}
		// Only the slack the last vector reads, not the whole buffer: it is as long as the longest subject so far.
		Arrays.fill(c, length, length + BLOCK, (byte) OUTSIDE);
		return c;
	}

	/**
	 * A bit per position of the block {@code c[block .. block + end)} whose cell has hits, and the words of the block
	 * in {@code bw}. {@link #LANES} positions at a time, without a branch on a presence bit.
	 */
	private static long blockMask(byte[] c, int block, int end, int[] present, int[] bw) {
		long mask = 0;
		for (int j = 0; j < end; j += LANES) {
			final int i = block + j;
			final IntVector c0 = widen(c, i), c1 = widen(c, i + 1), c2 = widen(c, i + 2), c3 = widen(c, i + 3),
					c4 = widen(c, i + 4);
			// One and tells whether any of the five letters is outside the alphabet: only those have that bit.
			final VectorMask<Integer> outside = c0.or(c1).or(c2).or(c3).or(c4).and(OUTSIDE)
					.compare(VectorOperators.NE, 0);
			// The cell of the word, its oldest letter lowest, as CompressedTable numbers them; 0 where a letter is
			// outside, so that the presence lookup below stays within present.
			final IntVector word = c0.add(c1.mul(15)).add(c2.mul(225)).add(c3.mul(3375)).add(c4.mul(50625))
					.blend(0, outside);
			word.intoArray(bw, j);
			// The presence bits one lane at a time, but still without a branch on them: a gather of the indices costs
			// more here than the loads it saves.
			long bits = 0;
			for (int k = 0; k < LANES; k++) {
				final int x = bw[j + k];
				bits |= (long) (present[x >>> 5] >>> x & 1) << k;
			}
			mask |= (bits & ~outside.toLong()) << j;
		}
		return mask;
	}

	/** {@link #LANES} compressed letters from {@code c[i]}, as ints. */
	private static IntVector widen(byte[] c, int i) {
		return (IntVector) ByteVector.fromArray(BYTES, c, i).convertShape(VectorOperators.B2I, INTS, 0);
	}

	/** {@link #positionsWithHits} as BLAST scans, one position at a time: for a native image. */
	private int scalarPositionsWithHits(byte[] residues, int start, int length, int[] present) {
		final byte[] compress = CompressedAlphabet.COMPRESS;
		final int[] w = words, pos = at;
		int n = 0;
		// The word ending at each position, its oldest letter lowest; good counts the letters since one outside the
		// alphabet.
		int word = 0, good = 0;
		for (int p = 0; p < length; p++) {
			final int c = compress[residues[start + p]];
			if (c >= CompressedAlphabet.SIZE) {
				good = 0;
				word = 0;
				continue;
			}
			word = word / CompressedAlphabet.SIZE + c * NEWEST;
			if (++good < WORD) {
				continue;
			}
			if ((present[word >>> 5] >>> word & 1) == 0) {
				continue; // no hits: cells not read
			}
			w[n] = word;
			pos[n++] = p - (WORD - 1);
		}
		return n;
	}
}
