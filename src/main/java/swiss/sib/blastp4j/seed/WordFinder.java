package swiss.sib.blastp4j.seed;

import java.util.List;

import swiss.sib.blastp4j.stats.Blosum62;

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
 * The scan makes no calls: the extensions are written out in it, and the alignments it finds go to a buffer that
 * {@link #find} turns into {@link UngappedHsp}s afterwards, in the order found.
 * <p>
 * Not thread safe: use one per thread.
 */
public final class WordFinder {
	/**
	 * BLOSUM62 with 32 columns per row, so that a row starts at a multiple of 32. Static, so that C2 knows its length and
	 * need not check the lookups against it.
	 */
	private static final byte[] PADDED = new byte[32 * 32];
	static {
		byte[] matrix = Blosum62.matrix();
		for (int a = 0; a < Blosum62.SIZE; a++) {
			for (int b = 0; b < Blosum62.SIZE; b++) {
				PADDED[a * 32 + b] = matrix[a * Blosum62.SIZE + b];
			}
		}
	}

	/** The ints per alignment in the buffer: query start, subject start, length, score, hit query and subject offsets. */
	private static final int FIELDS = 6;
	/** Set in a diagonal's state while hits on it are ignored, as they lie within the extension. */
	private static final int EXTENDED = 1 << 30;

	private final byte[] query;
	/** The start of each query residue's row of the padded matrix. */
	private final int[] rowOfQuery;
	private final LookupTable table;
	private final int window, xDrop, cutoff;
	/** The most hits one word can have, which bounds the alignments found in a subject. */
	private final int maxNeighbours;
	/**
	 * Per diagonal the subject offset of the last hit, with {@link #EXTENDED} set when the diagonal was extended to the
	 * right. One array, so that a hit reads and writes a single element.
	 */
	private final int[] diagState;
	private final int diagMask;
	private int diagOffset;
	/** The alignments found by the last scan, {@link #FIELDS} ints each. */
	private int[] found = new int[0];

	public WordFinder(byte[] query, LookupTable table, SeedParameters p) {
		this.query = query;
		this.table = table;
		this.window = p.window();
		this.xDrop = p.xDrop();
		this.cutoff = p.cutoff();
		int length = 1;
		while (length < query.length + window) {
			length <<= 1;
		}
		rowOfQuery = new int[query.length];
		for (int q = 0; q < query.length; q++) {
			rowOfQuery[q] = query[q] << 5;
		}
		int most = 0;
		for (int w = 0; w < LookupTable.WORDS; w++) {
			most = Math.max(most, table.start[w + 1] - table.start[w]);
		}
		maxNeighbours = most;
		diagState = new int[length];
		diagMask = length - 1;
		diagOffset = window;
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} and adds them to
	 * {@code out}, in the order BLAST finds them (BLAST later sorts them by score).
	 */
	public void find(byte[] residues, int start, int length, List<UngappedHsp> out) {
		final int[] cellStart = table.start, offsets = table.offsets;
		final int wordSize = LookupTable.WORD_SIZE;
		final int mask = LookupTable.WORDS - 1;
		final int[] diags = diagState;
		final int diagMaskL = diagMask, diagOffsetL = diagOffset, windowL = window;
		final int xDropL = xDrop, cutoffL = cutoff, queryLength = query.length;
		final int[] row = rowOfQuery;
		final byte[] matrix = PADDED;
		// Every hit of the subject can give one alignment at most.
		final int needed = FIELDS * length * maxNeighbours;
		if (found.length < needed) {
			found = new int[needed];
		}
		final int[] result = found;
		int n = 0;
		if (length >= wordSize) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - wordSize; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					final int q = offsets[k];
					final int diag = (q - s) & diagMaskL;
					final int state = diags[diag];
					if ((state & EXTENDED) == 0) {
						if (s - (state - diagOffsetL) >= windowL) {
							// The usual hit: the last hit on its diagonal is a window or more back.
							diags[diag] = s + diagOffsetL;
							continue;
						}
					} else {
						// Already extended past this hit: skip it; otherwise it starts a new pair.
						if (s + diagOffsetL >= (state & ~EXTENDED)) {
							diags[diag] = s + diagOffsetL;
						}
						continue;
					}
					final int last = state - diagOffsetL;
					if (s - last < wordSize) {
						continue; // overlaps the previous hit
					}

					// s_BlastAaExtendTwoHit: the second word ends the left extension's start; the right extension
					// runs only if the left one reached the first word (which ends at sLeft).
					final int sLeft = last + wordSize;
					int score = 0, bestWord = 0, rightD = 0;
					for (int i = 0; i < wordSize; i++) {
						score += matrix[(row[q + i] + residues[start + s + i]) & 1023];
						if (score > bestWord) {
							bestWord = score;
							rightD = i + 1;
						}
					}
					final int sRight = s + rightD, qRight = q + rightD;

					// s_BlastAaExtendLeft from just before (qRight, sRight), with a start score of 0.
					final int sOffL = sRight - 1, qOffL = qRight - 1;
					final int nL = Math.min(sOffL, qOffL);
					final int sBaseL = start + sOffL - nL, qBaseL = qOffL - nL;
					score = 0;
					int leftScore = 0, bestI = nL + 1, floorL = -xDropL;
					for (int i = nL; i >= 0; i--) {
						score += matrix[(row[qBaseL + i] + residues[sBaseL + i]) & 1023];
						if (score > leftScore) {
							leftScore = score;
							bestI = i;
							floorL = leftScore - xDropL;
						}
						if (score <= floorL) {
							break;
						}
					}
					final int leftD = nL - bestI + 1;

					// s_BlastAaExtendRight, when the left extension reached the first word.
					int rightScore = 0, rightLen = 0, lastSubjectOffset = sRight;
					final boolean rightExtended = leftD >= sRight - sLeft;
					if (rightExtended) {
						final int nR = Math.min(length - sRight, queryLength - qRight);
						int sc = leftScore, maxR = leftScore, bestJ = -1, floorR = Math.max(0, leftScore - xDropL);
						int i = 0;
						for (; i < nR; i++) {
							sc += matrix[(row[qRight + i] + residues[start + sRight + i]) & 1023];
							if (sc > maxR) {
								maxR = sc;
								bestJ = i;
								floorR = Math.max(0, maxR - xDropL);
							}
							if (sc <= floorR) {
								break;
							}
						}
						rightScore = maxR;
						rightLen = bestJ + 1;
						lastSubjectOffset = sRight + i;
					}

					final int best = Math.max(leftScore, rightScore);
					if (best >= cutoffL) {
						final int at = FIELDS * n++;
						result[at] = qRight - leftD;
						result[at + 1] = sRight - leftD;
						result[at + 2] = leftD + rightLen;
						result[at + 3] = best;
						result[at + 4] = q;
						result[at + 5] = s;
					}
					if (rightExtended) {
						diags[diag] = (lastSubjectOffset - (wordSize - 1) + diagOffsetL) | EXTENDED;
					} else {
						diags[diag] = s + diagOffsetL;
					}
				}
			}
		}
		for (int e = 0; e < n; e++) {
			final int at = FIELDS * e;
			out.add(new UngappedHsp(result[at], result[at + 1], result[at + 2], result[at + 3], result[at + 4],
					result[at + 5]));
		}
		// The next subject's positions start beyond everything stored so far.
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(diagState, 0);
		} else {
			diagOffset += length + window;
		}
	}
}
