package swiss.sib.blastp4j.seed;

import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.util.List;
import java.util.Map;

import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorShape;
import jdk.incubator.vector.VectorSpecies;
import swiss.sib.blastp4j.stats.Blosum62;
import swiss.sib.blastp4j.util.HiddenCopies;
import swiss.sib.blastp4j.util.LruCache;

/**
 * blastp's word finder ({@code s_BlastAaWordFinder_TwoHit} and {@code s_BlastAaExtendTwoHit} in {@code aa_ungapped.c}):
 * scans a subject for words in the lookup table and extends two hits on the same diagonal into an ungapped alignment.
 * <p>
 * Per diagonal it keeps the subject offset of the last hit. A hit is extended when the last hit on its diagonal is at
 * least a word length and less than the window behind it. The extension goes left from the end of the second word
 * with X-drop; only when it reaches the first word does it also go right. Once a diagonal has been extended to the
 * right, hits on it are ignored until they lie beyond the extension. Ungapped alignments scoring at least the cutoff
 * are reported.
 * <p>
 * The diagonal array has BLAST's size and indexing (a power of two of at least query length plus window, indexed by
 * query minus subject offset), so diagonals that share an entry behave as in BLAST. Instead of clearing it for each
 * subject, an offset is added to the stored positions, as BLAST does.
 * <p>
 * Not thread safe: use one per thread.
 */
public final class WordFinder implements WordScan {
	/**
	 * The copies by their constants, the most recently used {@code blastp4j.specializations} of them (256): the cutoff
	 * depends on the query and the database, so a long-running daemon would otherwise keep a copy for each one it
	 * met. A dropped copy is a hidden class no longer referenced, which the JVM can unload.
	 */
	private static final LruCache<Key, MethodHandle> COPIES = new LruCache<>(
			Integer.getInteger("blastp4j.specializations", 256));

	/** A copy's constants. */
	private record Key(int window, int xDrop, int cutoff, int diagMask) {
	}

	/**
	 * BLOSUM62 with 32 columns per row, so that a row starts at a multiple of 32. Static, so that C2 knows its length and
	 * need not check the lookups against it.
	 */
	static final byte[] PADDED = new byte[32 * 32];
	static {
		byte[] matrix = Blosum62.matrix();
		for (int a = 0; a < Blosum62.SIZE; a++) {
			for (int b = 0; b < Blosum62.SIZE; b++) {
				PADDED[a * 32 + b] = matrix[a * Blosum62.SIZE + b];
			}
		}
	}

	private final byte[] query;
	/** The start of each query residue's row of the padded matrix. */
	private final int[] rowOfQuery;
	private final LookupTable table;
	private final int window, xDrop, cutoff;
	/**
	 * Per diagonal the subject offset of the last hit, with {@link #EXTENDED} set when the diagonal was extended to the
	 * right. One array, so that a hit reads and writes a single element.
	 */
	private final int[] diagState;
	/** The offset of this query in the diagonals of {@link #diagState}: 0, or in a batch its offset in the concatenation. */
	private final int diagBase;
	/** Set in a diagonal's state while hits on it are ignored, as they lie within the extension. */
	static final int EXTENDED = 1 << 30;
	/**
	 * Whether the ungapped extensions track their best score without a branch on a new maximum (a select instead, which
	 * the predictor need not learn), so that only the X-drop exit branches: on AArch64 only. On the runners (654b5e5,
	 * run 42000813261) that took 0.8% less wall time on the Neoverse-N2, with 6% fewer branch misses, and 1.8% more on
	 * x86-64 (an EPYC 9V74).
	 */
	static final boolean BRANCH_FREE_MAX = "aarch64".equals(System.getProperty("os.arch"));
	private final int diagMask;
	private int diagOffset;
	// Output of the last extension.
	private int hspQuery, hspSubject, hspLength, lastSubjectOffset;
	private boolean rightExtended;

	/**
	 * The machine's preferred vector of ints: one subject position per lane. bin/jvm.options selects the width (AVX-512
	 * on Zen 5). Capped at 256 bits before: with AVX-512 (16 lanes) a single-query search of CDK1 took 10% longer on a
	 * Xeon 8573C (BlastpSearchBench, two runs), and the full search of 500 bean queries 8% longer on an EPYC 9V45.
	 */
	private static final VectorSpecies<Integer> INTS = IntVector.SPECIES_PREFERRED;
	/** Positions per vector. */
	private static final int LANES = INTS.length();
	/**
	 * The letters of {@link #LANES} positions, before they are widened to ints: the smallest byte vector with at least
	 * that many lanes (there is no byte vector smaller than 64 bits, so it can have more).
	 */
	private static final VectorSpecies<Byte> BYTES = VectorSpecies.of(byte.class,
			VectorShape.forBitSize(Math.max(64, LANES * Byte.SIZE)));
	/** Whether the byte vector has one lane per position, so that its loads stay within the subject. */
	private static final boolean BYTE_PER_LANE = BYTES.length() == LANES;
	/**
	 * Whether this is a GraalVM native image (the property is set there, also while WordFinder is initialized at build
	 * time). Graal does not compile a Vector API gather to a gather instruction, so there the table is read with scalar
	 * loads; C2 does, and is faster with it.
	 */
	private static final boolean NATIVE_IMAGE = System.getProperty("org.graalvm.nativeimage.imagecode") != null;
	/** Positions per block: one bit each in a long. */
	private static final int BLOCK = 64;
	/** The words of the current block of subject positions. */
	private final int[] blockWords = new int[BLOCK];
	/** The first two letters of the words of the current block, to index the presence table. */
	private final int[] blockPairs = new int[BLOCK];
	/**
	 * The lookup table coverage (expected fraction of subject positions with neighbours) from which {@link #find}, one
	 * position at a time, is faster than {@link #findBlocks}, a block of positions with vectors. The vectors pay by
	 * skipping positions without neighbours, so the fewer of those, the less they gain, and the narrower the vectors,
	 * the sooner they lose. With
	 * SeedSearchBench on GitHub's runners (run 41454777233), against the vector scan: with 8 lanes (x64, AVX2) EGFR
	 * (coverage 0.81) was 1.5% slower and BRCA1 (0.88) 2.0% faster one position at a time; with 4 lanes (AArch64, NEON)
	 * EGFR 5.7% and BRCA1 4.0% faster, FBF1 (0.71) the same.
	 */
	static final double VECTOR_MAX_COVERAGE = LANES >= 8 ? 0.85 : 0.75;
	/** Whether {@link #findBlocks} is for this query: its table has a coverage below {@link #VECTOR_MAX_COVERAGE}. */
	private final boolean vectorScan;

	public WordFinder(byte[] query, LookupTable table, SeedParameters p) {
		this(query, table, p, null, 0);
	}

	/**
	 * A word finder whose diagonals are those of the queries concatenated in {@code diags}, this one at offset
	 * {@code base} (for {@link BatchWordFinder}), or an array of its own when {@code diags} is null. Only {@link #hit}
	 * and {@link #endSubject} use the shared array: such a finder does not scan.
	 */
	WordFinder(byte[] query, LookupTable table, SeedParameters p, int[] diags, int base) {
		this.query = query;
		this.table = table;
		this.window = p.window();
		this.xDrop = p.xDrop();
		this.cutoff = p.cutoff();
		final int length = diags != null ? diags.length : diagonals(query.length, window);
		rowOfQuery = new int[query.length];
		for (int q = 0; q < query.length; q++) {
			rowOfQuery[q] = query[q] << 5;
		}
		diagState = diags != null ? diags : new int[length];
		diagBase = base;
		diagMask = length - 1;
		diagOffset = window;
		vectorScan = table.coverage() < VECTOR_MAX_COVERAGE;
	}

	/** The diagonal array's length for a query: a power of two of at least the query length plus the window. */
	static int diagonals(int queryLength, int window) {
		int length = 1;
		while (length < queryLength + window) {
			length <<= 1;
		}
		return length;
	}

	/**
	 * A word finder for this query whose window, X-drop, cutoff and diagonal mask are constants: a copy of this class,
	 * made with the class-file API and defined as a hidden class, with each read of those fields replaced by its
	 * value. One copy per combination of values, so queries that share them share a copy. In a native image, which
	 * cannot define classes at run time, or where the copy fails, a {@link WordFinder}.
	 */
	public static WordScan of(byte[] query, LookupTable table, SeedParameters p) {
		if (!NATIVE_IMAGE) {
			final Key key = new Key(p.window(), p.xDrop(), p.cutoff(), diagonals(query.length, p.window()) - 1);
			try {
				return (WordScan) COPIES.computeIfAbsent(key, k -> {
					try {
						return specialize(k);
					} catch (Throwable e) {
						throw new IllegalStateException(e);
					}
				}).invoke(query, table, p);
			} catch (Throwable e) {
				// Fall back to the fields.
			}
		}
		return new WordFinder(query, table, p);
	}

	/** This class with each read of the fields in {@code key} replaced by the value; its constructor. */
	private static MethodHandle specialize(Key key) throws Throwable {
		final MethodHandles.Lookup lookup = HiddenCopies.define(MethodHandles.lookup(), WordFinder.class,
				Map.of("window", key.window(), "xDrop", key.xDrop(), "cutoff", key.cutoff(), "diagMask",
						key.diagMask()));
		return lookup.findConstructor(lookup.lookupClass(),
				MethodType.methodType(void.class, byte[].class, LookupTable.class, SeedParameters.class));
	}

	/**
	 * Whether to call {@link #findBlocks} rather than {@link #find} for this query: its table has a coverage below
	 * {@link #VECTOR_MAX_COVERAGE}. The caller chooses, so that C2 inlines {@link #find} into the caller as before the
	 * vector scan: a find that chose between the scans was not inlined, and its scalar scan ran 5% more instructions.
	 */
	@Override
	public boolean vectorScan() {
		return vectorScan;
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} and adds them to
	 * {@code out}, in the order BLAST finds them (BLAST later sorts them by score). Scans one position at a time.
	 */
	@Override
	public void find(byte[] residues, int start, int length, List<UngappedHsp> out) {
		final int[] cellStart = table.start(), offsets = table.offsets();
		final int wordSize = LookupTable.WORD_SIZE;
		final int mask = LookupTable.WORDS - 1;
		// The diagonal fields as locals, for the usual hit below.
		final int[] diags = diagState;
		final int diagMaskL = diagMask, diagOffsetL = diagOffset, windowL = window;
		final int lastL = windowL - LookupTable.WORD_SIZE - 1;
		if (length >= wordSize) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - wordSize; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					final int q = offsets[k];
					final int diag = (q - s) & diagMaskL;
					final int state = diags[diag];
					// As in BatchWordFinder.find: stored or ignored here, with one branch, unless it is a hit to
					// extend or on an extended diagonal.
					final int diff = s - (state - diagOffsetL), x = diff - LookupTable.WORD_SIZE;
					if ((~(x | (lastL - x)) | state << 1) >= 0) {
						final int older = (diff - windowL) >> 31;
						diags[diag] = (state & older) | ((s + diagOffsetL) & ~older);
					} else {
						hit(residues, start, length, q, s, out);
					}
				}
			}
		}
		// The next subject's positions start beyond everything stored so far.
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(diagState, 0);
		} else {
			diagOffset += length + window;
		}
	}

	/** As {@link #find}, scanning blocks of positions and finding those whose words have neighbours with vectors. */
	@Override
	public void findBlocks(byte[] residues, int start, int length, List<UngappedHsp> out) {
		final int[] cellStart = table.start(), offsets = table.offsets(), presence = table.presence();
		final int[] words = blockWords, pairs = blockPairs;
		final int positions = length - LookupTable.WORD_SIZE + 1;
		// The diagonal fields as locals, for the usual hit below.
		final int[] diags = diagState;
		final int diagMaskL = diagMask, diagOffsetL = diagOffset, windowL = window;
		final int lastL = windowL - LookupTable.WORD_SIZE - 1;
		for (int block = 0; block < positions; block += BLOCK) {
			// First the words of the block and a bit per position whose word has neighbours: LANES positions at a
			// time with vectors, the rest one at a time, without branches.
			final int n = Math.min(BLOCK, positions - block);
			// A vector of positions from p on loads the bytes up to start + p + BYTES.length() + 1. With one byte per
			// lane those are letters of the subject; with more they go beyond it, and must still be in the array.
			final int vectorN = BYTE_PER_LANE ? n
					: Math.min(n, residues.length - start - BYTES.length() - 2 + LANES - block);
			long present = 0;
			int j = 0;
			for (; j + LANES <= vectorN; j += LANES) {
				// The letters of the LANES words, as ints. No vector helper method: when Graal does not inline one, it
				// runs the vector operations of the whole loop in Java (native image).
				final int i = start + block + j;
				IntVector a = (IntVector) ByteVector.fromArray(BYTES, residues, i)
						.convertShape(VectorOperators.B2I, INTS, 0);
				IntVector b = (IntVector) ByteVector.fromArray(BYTES, residues, i + 1)
						.convertShape(VectorOperators.B2I, INTS, 0);
				IntVector c = (IntVector) ByteVector.fromArray(BYTES, residues, i + 2)
						.convertShape(VectorOperators.B2I, INTS, 0);
				IntVector ab = a.lanewise(VectorOperators.LSHL, LookupTable.CHAR_BITS).or(b);
				ab.lanewise(VectorOperators.LSHL, LookupTable.CHAR_BITS).or(c).intoArray(words, j);
				if (!NATIVE_IMAGE) {
					// The presence bit of word abc is bit c of presence[ab]: one gather from a 4 KB table.
					ab.intoArray(pairs, j);
					present |= IntVector.fromArray(INTS, presence, 0, pairs, j).lanewise(VectorOperators.LSHR, c)
							.and(1).compare(VectorOperators.NE, 0).toLong() << j;
				}
			}
			if (NATIVE_IMAGE) {
				for (int k = 0; k < j; k++) {
					final int word = words[k];
					present |= (long) ((presence[word >>> LookupTable.CHAR_BITS] >>> word) & 1) << k;
				}
			}
			for (; j < n; j++) {
				final int word = LookupTable.word(residues, start + block + j);
				words[j] = word;
				present |= (long) ((presence[word >>> LookupTable.CHAR_BITS] >>> word) & 1) << j;
			}
			// Then the hits of those positions, in subject order.
			for (; present != 0; present &= present - 1) {
				final int b = Long.numberOfTrailingZeros(present);
				final int s = block + b, word = words[b];
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					final int q = offsets[k];
					final int diag = (q - s) & diagMaskL;
					final int state = diags[diag];
					// As in BatchWordFinder.find: stored or ignored here, with one branch, unless it is a hit to
					// extend or on an extended diagonal.
					final int diff = s - (state - diagOffsetL), x = diff - LookupTable.WORD_SIZE;
					if ((~(x | (lastL - x)) | state << 1) >= 0) {
						final int older = (diff - windowL) >> 31;
						diags[diag] = (state & older) | ((s + diagOffsetL) & ~older);
					} else {
						hit(residues, start, length, q, s, out);
					}
				}
			}
		}
		// The next subject's positions start beyond everything stored so far.
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(diagState, 0);
		} else {
			diagOffset += length + window;
		}
	}

	/** The end of a subject of {@code length} letters scanned elsewhere, as at the end of {@link #find}. */
	void endSubject(int length) {
		// The next subject's positions start beyond everything stored so far.
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(diagState, 0);
		} else {
			diagOffset += length + window;
		}
	}

	/**
	 * A word hit: query offset {@code q} at subject position {@code s}. Also called by {@link BatchWordFinder} for a
	 * hit that is not just stored.
	 */
	void hit(byte[] residues, int start, int length, int q, int s, List<UngappedHsp> out) {
		final int diag = (diagBase + q - s) & diagMask;
		final int state = diagState[diag];
		if ((state & EXTENDED) != 0) {
			// Already extended past this hit: skip it; otherwise it starts a new pair.
			if (s + diagOffset >= (state & ~EXTENDED)) {
				diagState[diag] = s + diagOffset;
			}
			return;
		}
		final int last = state - diagOffset;
		final int diff = s - last;
		if (diff >= window) {
			diagState[diag] = s + diagOffset;
			return;
		}
		if (diff < LookupTable.WORD_SIZE) {
			return; // overlaps the previous hit
		}
		if (q < diff) {
			// In a batch, the last hit on the diagonal was in a query before this one in the concatenation: as BLAST
			// (the check after BSearchContextInfo), this one starts a new pair. Never so for a query alone.
			diagState[diag] = s + diagOffset;
			return;
		}
		int score = extendTwoHit(residues, start, length, last + LookupTable.WORD_SIZE, s, q);
		if (score >= cutoff) {
			out.add(new UngappedHsp(hspQuery, hspSubject, hspLength, score, q, s));
		}
		if (rightExtended) {
			diagState[diag] = (lastSubjectOffset - (LookupTable.WORD_SIZE - 1) + diagOffset) | EXTENDED;
		} else {
			diagState[diag] = s + diagOffset;
		}
	}

	/**
	 * {@code s_BlastAaExtendTwoHit}: the second word ends the left extension's start; the right extension runs only
	 * if the left one reached the first word (which ends at {@code sLeft}).
	 */
	private int extendTwoHit(byte[] residues, int start, int length, int sLeft, int sRight, int qRight) {
		// One beyond the position within the word that gives the best starting score.
		int score = 0, leftScore = 0, rightD = 0;
		for (int i = 0; i < LookupTable.WORD_SIZE; i++) {
			score += PADDED[(rowOfQuery[qRight + i] + residues[start + sRight + i]) & 1023];
			if (score > leftScore) {
				leftScore = score;
				rightD = i + 1;
			}
		}
		qRight += rightD;
		sRight += rightD;
		rightExtended = false;
		lastSubjectOffset = sRight;

		// Extend left from just before (qRight, sRight).
		leftScore = extendLeft(residues, start, sRight - 1, qRight - 1);
		int leftD = extensionLength;
		int rightScore = 0;
		rightD = 0;
		if (leftD >= sRight - sLeft) {
			rightExtended = true;
			rightScore = extendRight(residues, start, length, sRight, qRight, leftScore);
			rightD = extensionLength;
		}
		hspQuery = qRight - leftD;
		hspSubject = sRight - leftD;
		hspLength = leftD + rightD;
		return Math.max(leftScore, rightScore);
	}

	private int extensionLength;

	/** {@code s_BlastAaExtendLeft} with a start score of 0. */
	private int extendLeft(byte[] residues, int start, int sOff, int qOff) {
		final int n = Math.min(sOff, qOff);
		int best = n + 1, score = 0, maxScore = 0, floor = -xDrop;
		final int sBase = start + sOff - n, qBase = qOff - n;
		for (int i = n; i >= 0; i--) {
			score += PADDED[(rowOfQuery[qBase + i] + residues[sBase + i]) & 1023];
			if (BRANCH_FREE_MAX) {
				// With masks, not score > maxScore ? i : best, which C2 compiles to a branch again (rarely taken).
				final int gt = (maxScore - score) >> 31;
				best = (i & gt) | (best & ~gt);
				maxScore = Math.max(maxScore, score);
				floor = maxScore - xDrop;
			} else if (score > maxScore) {
				maxScore = score;
				best = i;
				floor = maxScore - xDrop;
			}
			if (score <= floor) {
				break;
			}
		}
		extensionLength = n - best + 1;
		return maxScore;
	}

	/** {@code s_BlastAaExtendRight}, continuing from {@code maxScore}. */
	private int extendRight(byte[] residues, int start, int length, int sOff, int qOff, int maxScore) {
		final int n = Math.min(length - sOff, query.length - qOff);
		int best = -1, score = maxScore, i;
		int floor = Math.max(0, maxScore - xDrop);
		for (i = 0; i < n; i++) {
			score += PADDED[(rowOfQuery[qOff + i] + residues[start + sOff + i]) & 1023];
			if (BRANCH_FREE_MAX) {
				final int gt = (maxScore - score) >> 31; // as in extendLeft
				best = (i & gt) | (best & ~gt);
				maxScore = Math.max(maxScore, score);
				floor = Math.max(0, maxScore - xDrop);
			} else if (score > maxScore) {
				maxScore = score;
				best = i;
				floor = Math.max(0, maxScore - xDrop);
			}
			if (score <= floor) {
				break;
			}
		}
		extensionLength = best + 1;
		lastSubjectOffset = sOff + i;
		return maxScore;
	}
}
