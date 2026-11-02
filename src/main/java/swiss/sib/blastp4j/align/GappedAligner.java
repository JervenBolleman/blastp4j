package swiss.sib.blastp4j.align;

import java.util.Arrays;

import swiss.sib.blastp4j.stats.Blosum62;

/**
 * blastp's gapped X-drop extensions ({@code blast_gapalign.c}): from a start point, a gapped alignment is extended to
 * the left and to the right, each direction by dynamic programming over a band that grows and shrinks with the
 * X-drop criterion. Three variants:
 * <ul>
 * <li>{@link #semiGapped} ({@code Blast_SemiGappedAlign}): score only, used by the preliminary stage;</li>
 * <li>{@link #restricted} ({@code s_RestrictedGappedAlign}): score only, gaps may only start every 10 positions; an
 * approximation the preliminary stage uses for subjects whose best ungapped alignment is weak;</li>
 * <li>{@link #alignEx} ({@code ALIGN_EX}): with traceback, used by the final stage.</li>
 * </ul>
 * The code follows BLAST's closely, including its tie-breaking, so that the same alignments come out.
 * <p>
 * Sequences are NCBIstdaa codes. The subject is read from a database array, where it is surrounded by sentinels; the
 * dynamic programming reads one letter beyond each end of the subject, as BLAST does.
 * <p>
 * Not thread safe: use one per thread.
 */
public final class GappedAligner {
	static final int MININT = Integer.MIN_VALUE / 2;
	/** Gaps in the restricted alignment may only start at offsets that are multiples of this. */
	static final int RESTRICT_SIZE = 10;
	/** Traceback flags: the gap in A (or B) continues from the cell before. */
	static final int EXTEND_GAP_A = 0x10, EXTEND_GAP_B = 0x40;
	/** Window for the best start point of a gapped extension ({@code HSP_MAX_WINDOW}). */
	static final int HSP_MAX_WINDOW = 11;

	private static final int N = Blosum62.SIZE;
	/** The scores, {@code matrix[query * N + subject]}; BLOSUM62 unless set otherwise. */
	private final int[] matrix = new int[N * N];
	private final int gapOpen, gapExtend;

	private int[] best = new int[256], bestGap = new int[256];
	/** semiGapped's best and bestGap, interleaved. */
	private int[] dp = new int[512];
	/** semiGapped's B letters in the order it reads them. */
	private byte[] bLetters = new byte[256];
	private byte[] trace = new byte[1 << 16];
	private int[] rowStart = new int[256];

	/** Results of the last extension: offsets of the best cell (in letters from the start point). */
	int aOffset, bOffset;

	/** Results of {@link #extendScoreOnly} and {@link #extendWithTraceback}. */
	public int score, queryStart, queryStop, subjectStart, subjectStop;
	/** Result of {@link #extendWithTraceback}. */
	public EditScript editScript;

	private final EditScript left = new EditScript(), right = new EditScript();

	/**
	 * @param gapOpen   BLAST's gap open cost (11): a gap of k letters costs {@code gapOpen + k * gapExtend}
	 * @param gapExtend BLAST's gap extension cost (1)
	 */
	public GappedAligner(int gapOpen, int gapExtend) {
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		for (int a = 0; a < N; a++) {
			for (int b = 0; b < N; b++) {
				matrix[a * N + b] = Blosum62.score(a, b);
			}
		}
	}

	/** Uses {@code m[query][subject]} (28 x 28) as the score matrix from now on. */
	public void setMatrix(int[][] m) {
		for (int a = 0; a < N; a++) {
			System.arraycopy(m[a], 0, matrix, a * N, N);
		}
	}

	private void ensureDp(int n) {
		if (best.length < n) {
			int size = Math.max(n, best.length * 2);
			best = new int[size];
			bestGap = new int[size];
			dp = new int[2 * size];
		}
	}

	/**
	 * B's letters in the order semiGapped reads them: {@code b[bBase + 1 .. bBase + n + 1]} forwards,
	 * {@code b[bBase + n - 1 .. bBase - 1]} backwards (the last of each is a sentinel), cut at the ends of {@code b}.
	 */
	private byte[] bLetters(byte[] b, int bBase, int n, boolean reverse) {
		if (bLetters.length < n + 1) {
			bLetters = new byte[Math.max(n + 1, bLetters.length * 2)];
		}
		final byte[] letters = bLetters;
		if (!reverse) {
			System.arraycopy(b, bBase + 1, letters, 0, Math.min(n + 1, b.length - bBase - 1));
		} else {
			for (int k = 0, end = Math.min(n + 1, bBase + n); k < end; k++) {
				letters[k] = b[bBase + n - 1 - k];
			}
		}
		return letters;
	}

	/** Letter {@code a} (1-based) of A, read forwards from {@code aBase} or backwards from {@code aBase + m}. */
	private static int letterA(byte[] a, int aBase, int m, int index, boolean reverse) {
		return reverse ? a[aBase + m - index] : a[aBase + index];
	}

	/**
	 * Score-only X-drop extension ({@code Blast_SemiGappedAlign} with {@code score_only}). A is the query, B the
	 * subject. Forwards, A's letters are {@code a[aBase + 1 .. aBase + m]} and B's {@code b[bBase + 1 .. bBase + n]};
	 * backwards, they are {@code a[aBase + m - 1 .. aBase]} and {@code b[bBase + n - 1 .. bBase]}.
	 */
	int semiGapped(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop) {
		// The gap costs and the best cell in locals: the inner loop then reads and writes no fields.
		final int ext = gapExtend, goe = gapOpen + ext;
		xDrop = Math.max(xDrop, goe);
		aOffset = 0;
		bOffset = 0;
		if (n <= 0 || m <= 0) {
			return 0;
		}
		ensureDp(n + 2);
		// best and bestGap of a column next to each other, as BLAST's BlastGapDP: dp[2 * i] and dp[2 * i + 1].
		final int[] dp = this.dp;
		final int[] mat = matrix;
		// B read in one direction by one loop: its index is the loop's, so C2 checks its bounds once per row.
		final byte[] bs = bLetters(b, bBase, n, reverse);
		// The scores of A's current letter. Allocated here with a constant length, so that C2 knows an index masked
		// with 31 is in bounds and does not check it every cell.
		final int[] rowScores = new int[32];
		int score = -goe;
		dp[0] = 0;
		dp[1] = -goe;
		int i;
		for (i = 1; i <= n; i++) {
			if (score < -xDrop) {
				break;
			}
			dp[2 * i] = score;
			dp[2 * i + 1] = score - goe;
			score -= ext;
		}
		int bSize = i, bestScore = 0, first = 0, bestA = 0, bestB = 0;
		// bestScore - xDrop, kept up to date, so that the inner loop does not need xDrop.
		int floor = -xDrop;
		for (int ai = 1; ai <= m; ai++) {
			System.arraycopy(mat, letterA(a, aBase, m, ai, reverse) * N, rowScores, 0, N);
			score = MININT;
			int gapRow = MININT, last = first;
			for (int bi = first; bi < bSize; bi++) {
				int gapCol = dp[2 * bi + 1];
				int next = dp[2 * bi] + rowScores[bs[bi] & 31];
				score = Math.max(score, Math.max(gapCol, gapRow));
				if (score < floor) {
					if (bi == first) {
						first++;
					} else {
						dp[2 * bi] = MININT;
					}
				} else {
					last = bi;
					if (score > bestScore) {
						bestScore = score;
						floor = score - xDrop;
						bestA = ai;
						bestB = bi;
					}
					gapRow -= ext;
					gapCol -= ext;
					dp[2 * bi + 1] = Math.max(score - goe, gapCol);
					gapRow = Math.max(score - goe, gapRow);
					dp[2 * bi] = score;
				}
				score = next;
			}
			if (first == bSize) {
				break;
			}
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= floor && bSize <= n) {
					dp[2 * bSize] = gapRow;
					dp[2 * bSize + 1] = gapRow - goe;
					gapRow -= ext;
					bSize++;
				}
			}
			if (bSize <= n) {
				dp[2 * bSize] = MININT;
				dp[2 * bSize + 1] = MININT;
				bSize++;
			}
		}
		aOffset = bestA;
		bOffset = bestB;
		return bestScore;
	}

	/** {@code s_RestrictedGappedAlign}: as {@link #semiGapped}, but gaps only start every {@link #RESTRICT_SIZE}. */
	int restricted(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop) {
		final int goe = gapOpen + gapExtend;
		xDrop = Math.max(xDrop, goe);
		aOffset = 0;
		bOffset = 0;
		if (n <= 0 || m <= 0) {
			return 0;
		}
		ensureDp(n + 2);
		final int[] best = this.best, bestGap = this.bestGap;
		final int[] mat = matrix;
		int score = -goe;
		best[0] = 0;
		bestGap[0] = -goe;
		int i;
		for (i = 1; i <= n; i++) {
			if (score < -xDrop) {
				break;
			}
			best[i] = score;
			bestGap[i] = score - goe;
			score -= gapExtend;
		}
		int bSize = i, bestScore = 0, first = 0, bGap = 0;
		final int bStep = reverse ? -1 : 1;
		for (int ai = 1; ai <= m; ai++) {
			final int row = letterA(a, aBase, m, ai, reverse) * N;
			int bp = reverse ? bBase + n - first : bBase + first;
			score = MININT;
			int gapRow = MININT, last = first;
			if (ai % RESTRICT_SIZE != 0) {
				// A gap never starts in A in this row.
				for (int bi = first; bi < bSize; bi++) {
					bp += bStep;
					int next = best[bi] + mat[row + b[bp]];
					if (bi != bGap) {
						if (bestScore - score > xDrop) {
							best[bi] = MININT;
							if (bi == first) {
								first++;
							}
						} else {
							last = bi;
							if (score > bestScore) {
								bestScore = score;
								aOffset = ai;
								bOffset = bi;
							}
							best[bi] = score;
						}
					} else {
						bGap += RESTRICT_SIZE;
						int gapCol = bestGap[bi];
						if (score < gapCol) {
							score = gapCol;
						}
						if (bestScore - score > xDrop) {
							best[bi] = MININT;
							if (bi == first) {
								first++;
							}
						} else {
							last = bi;
							if (score > bestScore) {
								bestScore = score;
								aOffset = ai;
								bOffset = bi;
							}
							gapCol -= gapExtend;
							bestGap[bi] = Math.max(score - goe, gapCol);
							best[bi] = score;
						}
					}
					score = next;
				}
				gapRow = score;
			} else {
				for (int bi = first; bi < bSize; bi++) {
					bp += bStep;
					int next = best[bi] + mat[row + b[bp]];
					if (bi != bGap) {
						if (score < gapRow) {
							score = gapRow;
						}
						if (bestScore - score > xDrop) {
							best[bi] = MININT;
							if (bi == first) {
								first++;
							}
						} else {
							last = bi;
							if (score > bestScore) {
								bestScore = score;
								aOffset = ai;
								bOffset = bi;
							}
							gapRow -= gapExtend;
							gapRow = Math.max(score - goe, gapRow);
							best[bi] = score;
						}
					} else {
						bGap += RESTRICT_SIZE;
						int gapCol = bestGap[bi];
						if (score < gapCol) {
							score = gapCol;
						}
						if (score < gapRow) {
							score = gapRow;
						}
						if (bestScore - score > xDrop) {
							best[bi] = MININT;
							if (bi == first) {
								first++;
							}
						} else {
							last = bi;
							if (score > bestScore) {
								bestScore = score;
								aOffset = ai;
								bOffset = bi;
							}
							gapRow -= gapExtend;
							gapCol -= gapExtend;
							bestGap[bi] = Math.max(score - goe, gapCol);
							gapRow = Math.max(score - goe, gapRow);
							best[bi] = score;
						}
					}
					score = next;
				}
			}
			if (first == bSize) {
				break;
			}
			// The first offset of the next row where a gap in B may start.
			int r = first % RESTRICT_SIZE;
			bGap = first;
			if (r > 0) {
				bGap += RESTRICT_SIZE - r;
			}
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= bestScore - xDrop && bSize <= n) {
					best[bSize] = gapRow;
					bestGap[bSize] = gapRow - goe;
					gapRow -= gapExtend;
					bSize++;
				}
			}
			if (bSize <= n) {
				best[bSize] = MININT;
				bestGap[bSize] = MININT;
				bSize++;
			}
		}
		return bestScore;
	}

	private void ensureTrace(int n) {
		if (trace.length < n) {
			trace = Arrays.copyOf(trace, Math.max(n, trace.length * 2));
		}
	}

	/**
	 * {@code ALIGN_EX}: as {@link #semiGapped}, recording the traceback, then following it from the best cell back to
	 * the start point. The operations are added to {@code ops} in that order.
	 */
	int alignEx(byte[] a, int aBase, int m, byte[] b, int bBase, int n, boolean reverse, int xDrop, EditScript ops) {
		final int goe = gapOpen + gapExtend;
		xDrop = Math.max(xDrop, goe);
		aOffset = 0;
		bOffset = 0;
		if (n <= 0 || m <= 0) {
			return 0;
		}
		ensureDp(n + 2);
		if (rowStart.length < m + 1) {
			rowStart = new int[Math.max(m + 1, rowStart.length * 2)];
		}
		final int[] best = this.best, bestGap = this.bestGap;
		final int[] mat = matrix;
		final int numExtraCells = gapExtend > 0 ? xDrop / gapExtend + 3 : n + 3;

		ensureTrace(n + 3);
		byte[] tb = trace;
		rowStart[0] = 0;
		int score = -goe;
		best[0] = 0;
		bestGap[0] = -goe;
		int i;
		for (i = 1; i <= n; i++) {
			if (score < -xDrop) {
				break;
			}
			best[i] = score;
			bestGap[i] = score - goe;
			score -= gapExtend;
			tb[i] = EditScript.DEL;
		}
		int used = i + 1;
		int bSize = i, bestScore = 0, first = 0;
		final int bStep = reverse ? -1 : 1;
		for (int ai = 1; ai <= m; ai++) {
			// Room for this row: from 'first' to at most numExtraCells beyond the previous row's end.
			ensureTrace(used + 2 + (bSize - first) + numExtraCells + 1);
			tb = trace;
			final int base = used + 1 - first; // tb[base + bi] is cell (ai, bi)
			rowStart[ai] = base;
			final int origFirst = first;
			final int row = letterA(a, aBase, m, ai, reverse) * N;
			int bp = reverse ? bBase + n - first : bBase + first;
			score = MININT;
			int gapRow = MININT, last = first;
			int bi;
			for (bi = first; bi < bSize; bi++) {
				bp += bStep;
				int gapCol = bestGap[bi];
				int next = best[bi] + mat[row + b[bp]];
				// The best of the three moves, and which one (ties keep the earlier one), without branches on the scores.
				int script = gapCol > score ? EditScript.INS : EditScript.SUB;
				final int afterCol = Math.max(score, gapCol);
				script = gapRow > afterCol ? EditScript.DEL : script;
				score = Math.max(afterCol, gapRow);
				if (bestScore - score > xDrop) {
					if (first == bi) {
						first++;
					} else {
						best[bi] = MININT;
					}
				} else {
					last = bi;
					if (score > bestScore) {
						bestScore = score;
						aOffset = ai;
						bOffset = bi;
					}
					gapRow -= gapExtend;
					gapCol -= gapExtend;
					if (gapCol < score - goe) {
						bestGap[bi] = score - goe;
					} else {
						bestGap[bi] = gapCol;
						script += EXTEND_GAP_B;
					}
					if (gapRow < score - goe) {
						gapRow = score - goe;
					} else {
						script += EXTEND_GAP_A;
					}
					best[bi] = score;
				}
				score = next;
				tb[base + bi] = (byte) script;
			}
			if (first == bSize) {
				break;
			}
			if (last < bSize - 1) {
				bSize = last + 1;
			} else {
				while (gapRow >= bestScore - xDrop && bSize <= n) {
					best[bSize] = gapRow;
					bestGap[bSize] = gapRow - goe;
					gapRow -= gapExtend;
					if (base + bSize >= tb.length) {
						ensureTrace(base + bSize + numExtraCells + 2);
						tb = trace;
					}
					tb[base + bSize] = EditScript.DEL;
					bSize++;
				}
			}
			used += Math.max(bi, bSize) - origFirst + 1;
			if (bSize <= n) {
				best[bSize] = MININT;
				bestGap[bSize] = MININT;
				bSize++;
			}
		}

		// Follow the traceback from the best cell to the start.
		int ai = aOffset, bi = bOffset;
		int script = EditScript.SUB;
		while (ai > 0 || bi > 0) {
			int next = tb[rowStart[ai] + bi];
			switch (script) {
			case EditScript.DEL:
				script = next & 0x07;
				if ((next & EXTEND_GAP_A) != 0) {
					script = EditScript.DEL;
				}
				break;
			case EditScript.INS:
				script = next & 0x07;
				if ((next & EXTEND_GAP_B) != 0) {
					script = EditScript.INS;
				}
				break;
			default:
				script = next & 0x07;
				break;
			}
			if (script == EditScript.DEL) {
				bi--;
			} else if (script == EditScript.INS) {
				ai--;
			} else {
				ai--;
				bi--;
			}
			ops.add(script, 1);
		}
		return bestScore;
	}

	/**
	 * {@code s_BlastProtGappedAlignment}: score-only extension of the query against the subject
	 * {@code residues[subjectBase .. subjectBase + subjectLength)} from the start point ({@code qOff}, {@code sOff}),
	 * which belongs to the left extension. Sets {@link #score} and the bounds (stops are exclusive).
	 */
	public void extendScoreOnly(byte[] query, byte[] residues, int subjectBase, int subjectLength, int qOff, int sOff,
			boolean restricted, int xDrop) {
		final int qLength = qOff + 1, sLength = sOff + 1;
		final int queryLength = query.length;
		int scoreLeft = 0, scoreRight = 0;
		// Left, including the start point.
		scoreLeft = restricted
				? restricted(query, 0, qLength, residues, subjectBase, sLength, true, xDrop)
				: semiGapped(query, 0, qLength, residues, subjectBase, sLength, true, xDrop);
		queryStart = qLength - aOffset;
		subjectStart = sLength - bOffset;
		if (qLength < queryLength && sLength < subjectLength) {
			scoreRight = restricted
					? restricted(query, qOff, queryLength - qLength, residues, subjectBase + sOff,
							subjectLength - sLength, false, xDrop)
					: semiGapped(query, qOff, queryLength - qLength, residues, subjectBase + sOff,
							subjectLength - sLength, false, xDrop);
			queryStop = aOffset + qOff + 1;
			subjectStop = bOffset + sOff + 1;
		} else {
			queryStop = qOff;
			subjectStop = sOff;
		}
		score = scoreLeft + scoreRight;
	}

	/**
	 * {@code BLAST_GappedAlignmentWithTraceback}: the final extension from ({@code qStart}, {@code sStart}), with
	 * traceback. Leading and trailing gaps of the resulting edit script are removed, with their cost added back, as
	 * BLAST does.
	 */
	public void extendWithTraceback(byte[] query, byte[] residues, int subjectBase, int subjectLength, int qStart,
			int sStart, int xDrop) {
		final int queryLength = query.length;
		left.clear();
		right.clear();
		int scoreLeft = alignEx(query, 0, qStart + 1, residues, subjectBase, sStart + 1, true, xDrop, left);
		queryStart = qStart - aOffset + 1;
		subjectStart = sStart - bOffset + 1;
		int scoreRight = 0;
		if (qStart < queryLength && sStart < subjectLength) {
			scoreRight = alignEx(query, qStart, queryLength - qStart - 1, residues, subjectBase + sStart,
					subjectLength - sStart - 1, false, xDrop, right);
			queryStop = qStart + aOffset + 1;
			subjectStop = sStart + bOffset + 1;
		} else {
			queryStop = qStart - 1;
			subjectStop = sStart - 1;
		}
		EditScript e = EditScript.join(left, right);
		while (e.size > 0 && e.ops[0] != EditScript.SUB) {
			scoreLeft += gapOpen + e.nums[0] * gapExtend;
			if (e.ops[0] == EditScript.DEL) {
				subjectStart += e.nums[0];
			} else {
				queryStart += e.nums[0];
			}
			e.removeFirst();
		}
		while (e.size > 0 && e.ops[e.size - 1] != EditScript.SUB) {
			scoreRight += gapOpen + e.nums[e.size - 1] * gapExtend;
			if (e.ops[e.size - 1] == EditScript.DEL) {
				subjectStop -= e.nums[e.size - 1];
			} else {
				queryStop -= e.nums[e.size - 1];
			}
			e.size--;
		}
		editScript = e;
		score = scoreLeft + scoreRight;
	}

	/**
	 * {@code BlastGetStartForGappedAlignment}: the query offset where the gapped extension of an ungapped alignment
	 * starts, the middle of its best-scoring window of 11 letters (or its middle, when it is that short).
	 */
	public int startForGapped(byte[] query, byte[] residues, int subjectBase, int qStart, int qLength, int sStart,
			int sLength) {
		if (qLength <= HSP_MAX_WINDOW) {
			return qStart + qLength / 2;
		}
		final int[] mat = matrix;
		int hspEnd = qStart + HSP_MAX_WINDOW;
		int q = qStart, s = subjectBase + sStart;
		int score = 0;
		for (int i = qStart; i < hspEnd; i++) {
			score += mat[query[q] * N + residues[s]];
			q++;
			s++;
		}
		int maxScore = score, maxOffset = hspEnd - 1;
		hspEnd = qStart + Math.min(qLength, sLength);
		for (int i = qStart + HSP_MAX_WINDOW; i < hspEnd; i++) {
			score -= mat[query[q - HSP_MAX_WINDOW] * N + residues[s - HSP_MAX_WINDOW]];
			score += mat[query[q] * N + residues[s]];
			if (score > maxScore) {
				maxScore = score;
				maxOffset = i;
			}
			q++;
			s++;
		}
		return maxScore > 0 ? maxOffset - HSP_MAX_WINDOW / 2 : qStart;
	}

	/**
	 * {@code BlastGetOffsetsForGappedAlignment}: the start point for the traceback of an HSP whose recorded start is
	 * (0, 0). Returns {@code {q, s}}, or null when no window scores above zero.
	 */
	public int[] offsetsForGapped(byte[] query, byte[] residues, int subjectBase, int qOffset, int qEnd, int sOffset,
			int sEnd) {
		final int[] mat = matrix;
		int qLength = qEnd - qOffset, sLength = sEnd - sOffset;
		if (qLength <= HSP_MAX_WINDOW) {
			return new int[] { qOffset + qLength / 2, sOffset + qLength / 2 };
		}
		int hspEnd = qOffset + HSP_MAX_WINDOW;
		int q = qOffset, s = subjectBase + sOffset;
		int score = 0;
		for (int i = qOffset; i < hspEnd; i++) {
			score += mat[query[q] * N + residues[s]];
			q++;
			s++;
		}
		int maxScore = score, maxOffset = hspEnd - 1;
		hspEnd = qOffset + Math.min(qLength, sLength);
		for (int i = qOffset + HSP_MAX_WINDOW; i < hspEnd; i++) {
			score -= mat[query[q - HSP_MAX_WINDOW] * N + residues[s - HSP_MAX_WINDOW]];
			score += mat[query[q] * N + residues[s]];
			if (score > maxScore) {
				maxScore = score;
				maxOffset = i;
			}
			q++;
			s++;
		}
		if (maxScore > 0) {
			return new int[] { maxOffset, maxOffset - qOffset + sOffset };
		}
		score = 0;
		q = qOffset + qLength - HSP_MAX_WINDOW;
		s = subjectBase + sOffset + sLength - HSP_MAX_WINDOW;
		for (int i = qEnd - HSP_MAX_WINDOW; i < qEnd; i++) {
			score += mat[query[q] * N + residues[s]];
			q++;
			s++;
		}
		if (score > 0) {
			return new int[] { qEnd - HSP_MAX_WINDOW / 2, sEnd - HSP_MAX_WINDOW / 2 };
		}
		return null;
	}
}
