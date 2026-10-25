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
 * Not thread safe: use one per thread.
 */
public final class WordFinder {
	private final byte[] query;
	private final byte[] matrix = Blosum62.matrix();
	private final LookupTable table;
	private final int window, xDrop, cutoff;
	private final int[] lastHit;
	private final boolean[] extended;
	private final int diagMask;
	private int diagOffset;
	// Output of the last extension.
	private int hspQuery, hspSubject, hspLength, lastSubjectOffset;
	private boolean rightExtended;

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
		lastHit = new int[length];
		extended = new boolean[length];
		diagMask = length - 1;
		diagOffset = window;
	}

	/**
	 * Finds the ungapped alignments of the subject {@code residues[start .. start + length)} and adds them to
	 * {@code out}, in the order BLAST finds them (BLAST later sorts them by score).
	 */
	public void find(byte[] residues, int start, int length, List<UngappedHsp> out) {
		final int[] cellStart = table.start, offsets = table.offsets;
		final long[] occupied = table.occupied;
		final int wordSize = LookupTable.WORD_SIZE;
		final int mask = LookupTable.WORDS - 1;
		if (length >= wordSize) {
			int word = (residues[start] << LookupTable.CHAR_BITS) | residues[start + 1];
			for (int s = 0; s <= length - wordSize; s++) {
				word = ((word << LookupTable.CHAR_BITS) | residues[start + s + 2]) & mask;
				if ((occupied[word >>> 6] & (1L << word)) == 0) {
					continue;
				}
				for (int k = cellStart[word], end = cellStart[word + 1]; k < end; k++) {
					hit(residues, start, length, offsets[k], s, out);
				}
			}
		}
		// The next subject's positions start beyond everything stored so far.
		if (diagOffset >= Integer.MAX_VALUE / 4) {
			diagOffset = window;
			java.util.Arrays.fill(lastHit, -window);
			java.util.Arrays.fill(extended, false);
		} else {
			diagOffset += length + window;
		}
	}

	private void hit(byte[] residues, int start, int length, int q, int s, List<UngappedHsp> out) {
		final int diag = (q - s) & diagMask;
		if (extended[diag]) {
			// Already extended past this hit: skip it; otherwise it starts a new pair.
			if (s + diagOffset >= lastHit[diag]) {
				lastHit[diag] = s + diagOffset;
				extended[diag] = false;
			}
			return;
		}
		final int last = lastHit[diag] - diagOffset;
		final int diff = s - last;
		if (diff >= window) {
			lastHit[diag] = s + diagOffset;
			return;
		}
		if (diff < LookupTable.WORD_SIZE) {
			return; // overlaps the previous hit
		}
		int score = extendTwoHit(residues, start, length, last + LookupTable.WORD_SIZE, s, q);
		if (score >= cutoff) {
			out.add(new UngappedHsp(hspQuery, hspSubject, hspLength, score, q, s));
		}
		if (rightExtended) {
			extended[diag] = true;
			lastHit[diag] = lastSubjectOffset - (LookupTable.WORD_SIZE - 1) + diagOffset;
		} else {
			lastHit[diag] = s + diagOffset;
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
			score += matrix[query[qRight + i] * Blosum62.SIZE + residues[start + sRight + i]];
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
		final byte[] m = matrix;
		final int n = Math.min(sOff, qOff);
		int best = n + 1, score = 0, maxScore = 0;
		final int sBase = start + sOff - n, qBase = qOff - n;
		for (int i = n; i >= 0; i--) {
			score += m[query[qBase + i] * Blosum62.SIZE + residues[sBase + i]];
			if (score > maxScore) {
				maxScore = score;
				best = i;
			}
			if (maxScore - score >= xDrop) {
				break;
			}
		}
		extensionLength = n - best + 1;
		return maxScore;
	}

	/** {@code s_BlastAaExtendRight}, continuing from {@code maxScore}. */
	private int extendRight(byte[] residues, int start, int length, int sOff, int qOff, int maxScore) {
		final byte[] m = matrix;
		final int n = Math.min(length - sOff, query.length - qOff);
		int best = -1, score = maxScore, i;
		for (i = 0; i < n; i++) {
			score += m[query[qOff + i] * Blosum62.SIZE + residues[start + sOff + i]];
			if (score > maxScore) {
				maxScore = score;
				best = i;
			}
			if (score <= 0 || maxScore - score >= xDrop) {
				break;
			}
		}
		extensionLength = best + 1;
		lastSubjectOffset = sOff + i;
		return maxScore;
	}
}
