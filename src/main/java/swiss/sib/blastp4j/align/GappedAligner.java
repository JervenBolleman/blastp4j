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
	/** The score matrix; BLOSUM62 unless set otherwise. */
	private ScoreMatrix matrix = ScoreMatrix.BLOSUM62;
	private final int gapOpen, gapExtend;

	private int[] best = new int[256], bestGap = new int[256];
	/** ExRows' copy of a band of {@link #best} and {@link #bestGap} before rows computed without traceback. */
	private int[] savedBest = new int[256], savedGap = new int[256];
	/** semiGapped's best and bestGap, interleaved. */
	private int[] dp = new int[512];
	/** semiGapped's B letters in the order it reads them. */
	private byte[] bLetters = new byte[256];
	private byte[] trace = new byte[1 << 16];
	private int[] rowStart = new int[256];

	/** Results of the last extension: offsets of the best cell (in letters from the start point). */
	int aOffset, bOffset;
	/** Results of {@link ExRows}' last row: its first live cell, its gap in A at the end, and the best score so far. */
	private int rowFirst, rowGapRow, rowBest;

	/** Results of {@link #extendScoreOnly} and {@link #extendWithTraceback}. */
	public int score, queryStart, queryStop, subjectStart, subjectStop;
	/** Result of {@link #extendWithTraceback}. */
	public EditScript editScript;

	private final EditScript left = new EditScript(), right = new EditScript();
	/**
	 * The row loops, with this aligner's gap costs as constants where that can be done ({@link Rows#of}): with BLOSUM62
	 * as constants too, and reading the matrix field for any other matrix.
	 */
	private final RowLoops rowLoops, rowLoopsBlosum62;
	/** The row loop for {@link #matrix}. */
	private RowLoops rows;
	/** alignEx's row loop, with this aligner's gap costs as constants where that can be done ({@link ExRows#of}). */
	private final ExRowLoops exRows;

	/**
	 * @param gapOpen   BLAST's gap open cost (11): a gap of k letters costs {@code gapOpen + k * gapExtend}
	 * @param gapExtend BLAST's gap extension cost (1)
	 */
	public GappedAligner(int gapOpen, int gapExtend) {
		this(gapOpen, gapExtend, Rows.of(gapOpen + gapExtend, gapExtend, null),
				Rows.of(gapOpen + gapExtend, gapExtend, Rows.BLOSUM62));
	}

	/** With the given row loops (the tests compare the copies with {@link Rows} itself). */
	GappedAligner(int gapOpen, int gapExtend, RowLoops rowLoops, RowLoops rowLoopsBlosum62) {
		this(gapOpen, gapExtend, rowLoops, rowLoopsBlosum62, ExRows.of(gapOpen + gapExtend, gapExtend));
	}

	/** With the given row loops, alignEx's too (the tests compare the copies with {@link ExRows} itself). */
	GappedAligner(int gapOpen, int gapExtend, RowLoops rowLoops, RowLoops rowLoopsBlosum62, ExRowLoops exRows) {
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.rowLoops = rowLoops;
		this.rowLoopsBlosum62 = rowLoopsBlosum62;
		this.rows = rowLoopsBlosum62;
		this.exRows = exRows;
	}

	public int gapOpen() {
		return gapOpen;
	}

	public int gapExtend() {
		return gapExtend;
	}

	/** Uses {@code m} as the score matrix from now on. */
	public void setMatrix(ScoreMatrix m) {
		matrix = m;
		rows = m == ScoreMatrix.BLOSUM62 ? rowLoopsBlosum62 : rowLoops;
	}

	private void ensureDp(int n) {
		if (best.length < n) {
			int size = Math.max(n, best.length * 2);
			best = new int[size];
			bestGap = new int[size];
			savedBest = new int[size];
			savedGap = new int[size];
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
		final int goe = gapOpen + gapExtend;
		xDrop = Math.max(xDrop, goe);
		aOffset = 0;
		bOffset = 0;
		if (n <= 0 || m <= 0) {
			return 0;
		}
		ensureDp(n + 2);
		// best and bestGap of a column next to each other, as BLAST's BlastGapDP: dp[2 * i] and dp[2 * i + 1].
		final int[] dp = this.dp;
		// B read in one direction by one loop: its index is the loop's, so C2 checks its bounds once per row.
		bLetters(b, bBase, n, reverse);
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
			score -= gapExtend;
		}
		return rows.rows(this, a, aBase, m, n, reverse, xDrop, i);
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
		final int[] mat = matrix.scores;
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
		// B read in one direction by one loop: its index is the loop's, so C2 checks its bounds once per row.
		final byte[] bs = bLetters(b, bBase, n, reverse);

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
		// The rest in a copy of ExRows with this aligner's gap costs as constants (ExRows.of).
		return exRows.rows(this, a, aBase, m, bs, n, reverse, xDrop, i, ops);
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
		trimEndGaps(left, scoreLeft, right, scoreRight);
	}

	/**
	 * The end of {@link #extendWithTraceback}, from its two alignEx results (the right one 0, 0, 0 and no operations
	 * where it does not run): sets {@link #score}, the bounds and {@link #editScript}. Also for alignEx's results
	 * computed elsewhere ({@link TracebackLanes}).
	 */
	public void joinTraceback(int queryLength, int subjectLength, int qStart, int sStart, int scoreLeft, int aLeft,
			int bLeft, EditScript left, int scoreRight, int aRight, int bRight, EditScript right) {
		queryStart = qStart - aLeft + 1;
		subjectStart = sStart - bLeft + 1;
		if (qStart < queryLength && sStart < subjectLength) {
			queryStop = qStart + aRight + 1;
			subjectStop = sStart + bRight + 1;
		} else {
			queryStop = qStart - 1;
			subjectStop = sStart - 1;
		}
		trimEndGaps(left, scoreLeft, right, scoreRight);
	}

	/**
	 * Joins the two sides' edit scripts into {@link #editScript} and sets {@link #score}, with end gaps trimmed (their
	 * cost added back) as BLAST does.
	 */
	private void trimEndGaps(EditScript left, int scoreLeft, EditScript right, int scoreRight) {
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
		final int[] mat = matrix.scores;
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
		final int[] mat = matrix.scores;
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

	/** A row loop copy's constants: the gap costs and the matrix (by identity; only {@link Rows#BLOSUM62} for now). */
	record Key(int goe, int ext, int[] matrix) {
		@Override
		public boolean equals(Object o) {
			return o instanceof Key k && k.goe == goe && k.ext == ext && k.matrix == matrix;
		}

		@Override
		public int hashCode() {
			return 31 * (31 * goe + ext) + System.identityHashCode(matrix);
		}
	}

	/**
	 * A copy of the row loop class {@code c} ({@link Rows} or {@link ExRows}), made with the class-file API and defined
	 * as a hidden class in this nest (so it reads the aligner's private fields as {@code c} does), with each read of
	 * {@code goe} and {@code ext} replaced by the constant and, with a matrix, each read of {@code MATRIX} by the hidden
	 * class's class data, a dynamic constant. Returns an instance, made with the constructor of the two gap costs.
	 */
	private static Object copy(Class<?> c, int goe, int ext, int[] matrix) throws Throwable {
		final byte[] bytes;
		try (java.io.InputStream in = c.getResourceAsStream(c.getName().substring(c.getName().lastIndexOf('.') + 1)
				+ ".class")) {
			bytes = in.readAllBytes();
		}
		final java.lang.classfile.ClassFile cf = java.lang.classfile.ClassFile.of();
		final java.lang.classfile.ClassModel model = cf.parse(bytes);
		final java.lang.constant.ClassDesc self = model.thisClass().asSymbol();
		final java.lang.classfile.CodeTransform constants = (code, element) -> {
			if (element instanceof java.lang.classfile.instruction.FieldInstruction f
					&& f.opcode() == java.lang.classfile.Opcode.GETFIELD && f.owner().asSymbol().equals(self)
					&& (f.name().equalsString("goe") || f.name().equalsString("ext"))) {
				code.pop(); // the receiver
				code.loadConstant(f.name().equalsString("goe") ? goe : ext);
			} else if (matrix != null && element instanceof java.lang.classfile.instruction.FieldInstruction f
					&& f.opcode() == java.lang.classfile.Opcode.GETSTATIC && f.owner().asSymbol().equals(self)
					&& f.name().equalsString("MATRIX")) {
				code.loadConstant(java.lang.constant.DynamicConstantDesc.ofNamed(
						java.lang.constant.ConstantDescs.BSM_CLASS_DATA, java.lang.constant.ConstantDescs.DEFAULT_NAME,
						java.lang.constant.ConstantDescs.CD_int.arrayType()));
			} else {
				code.with(element);
			}
		};
		// A hidden class takes its nest from the lookup, not from these attributes.
		final byte[] copy = cf.transformClass(model, java.lang.classfile.ClassTransform
				.dropping(e -> e instanceof java.lang.classfile.attribute.NestHostAttribute
						|| e instanceof java.lang.classfile.attribute.InnerClassesAttribute)
				.andThen(java.lang.classfile.ClassTransform.transformingMethodBodies(constants)));
		final java.lang.invoke.MethodHandles.Lookup lookup = matrix == null
				? java.lang.invoke.MethodHandles.lookup().defineHiddenClass(copy, true,
						java.lang.invoke.MethodHandles.Lookup.ClassOption.NESTMATE)
				: java.lang.invoke.MethodHandles.lookup().defineHiddenClassWithClassData(copy, matrix.clone(), true,
						java.lang.invoke.MethodHandles.Lookup.ClassOption.NESTMATE);
		return lookup.findConstructor(lookup.lookupClass(),
				java.lang.invoke.MethodType.methodType(void.class, int.class, int.class)).invoke(goe, ext);
	}

	/** The row loop of {@link #alignEx}: {@link ExRows}, or a copy of it with constant gap costs. */
	interface ExRowLoops {
		int rows(GappedAligner g, byte[] a, int aBase, int m, byte[] bs, int n, boolean reverse, int xDrop, int bSize,
				EditScript ops);
	}

	/**
	 * alignEx's rows of A after the initial row, and the traceback, with the gap costs as fields: {@code goe} (gapOpen +
	 * gapExtend) and {@code ext} (gapExtend), read once into locals. As for {@link Rows}, {@link #of} gives a copy in
	 * which they are constants, so that C2 folds them into the loop: BLAST's (11, 1), and with composition-based
	 * statistics, the default, the same scaled by 32.
	 */
	static final class ExRows implements ExRowLoops {
		/**
		 * Rows with traceback after the best cell's row, before the rows that follow are computed without it, score only
		 * (see {@link #rows}).
		 */
		static final int AHEAD = 2;

		/**
		 * Rows {@code from} onwards as {@link #rows} computes them, with the band {@code [first, bSize)} of the row before,
		 * scores only: the row of the first cell above {@code bestScore}, or -1 where the band ends (or A) before one.
		 * Changes {@code best} and {@code bestGap} (the caller keeps a copy), nothing else.
		 */
		private int scoreRows(GappedAligner g, byte[] a, int aBase, int m, byte[] bs, int n, boolean reverse, int xDrop,
				int from, int first, int bSize, int bestScore, int[] mat, int mask) {
			final int goe = this.goe, ext = this.ext;
			final int[] best = g.best, bestGap = g.bestGap;
			final int floor = bestScore - xDrop;
			for (int ai = from; ai <= m; ai++) {
				final int row = letterA(a, aBase, m, ai, reverse) << 5;
				int score = MININT, gapRow = MININT, last = first;
				for (int bi = first; bi < bSize; bi++) {
					int gapCol = bestGap[bi];
					final int next = best[bi] + mat[(row + bs[bi]) & mask];
					score = Math.max(Math.max(score, gapCol), gapRow);
					if (score < floor) {
						if (first == bi) {
							first++;
						} else {
							best[bi] = MININT;
						}
					} else {
						if (score > bestScore) {
							return ai;
						}
						last = bi;
						gapRow -= ext;
						gapCol -= ext;
						bestGap[bi] = Math.max(score - goe, gapCol);
						gapRow = Math.max(score - goe, gapRow);
						best[bi] = score;
					}
					score = next;
				}
				if (first == bSize) {
					return -1;
				}
				if (last < bSize - 1) {
					bSize = last + 1;
				} else {
					while (gapRow >= floor && bSize <= n) {
						best[bSize] = gapRow;
						bestGap[bSize] = gapRow - goe;
						gapRow -= ext;
						bSize++;
					}
				}
				if (bSize <= n) {
					best[bSize] = MININT;
					bestGap[bSize] = MININT;
					bSize++;
				}
			}
			return -1;
		}

		private static final java.util.Map<Key, ExRowLoops> COPIES = new java.util.concurrent.ConcurrentHashMap<>();

		private final int goe, ext;

		ExRows(int goe, int ext) {
			this.goe = goe;
			this.ext = ext;
		}

		/** The row loop for these gap costs: a copy with them as constants, or this class where that fails. */
		static ExRowLoops of(int goe, int ext) {
			if (Rows.NATIVE_IMAGE) {
				return new ExRows(goe, ext);
			}
			return COPIES.computeIfAbsent(new Key(goe, ext, null), key -> {
				try {
					return (ExRowLoops) copy(ExRows.class, goe, ext, null);
				} catch (Throwable e) {
					return new ExRows(goe, ext);
				}
			});
		}

		/**
		 * The cells {@code [first, bSize)} of row {@code ai} of {@link #rows}, with their traceback bytes at
		 * {@code tb[base + bi]}: its row's last live cell; the row's first live cell, its gap in A at the end and the best
		 * score in {@code g.rowFirst}, {@code g.rowGapRow} and {@code g.rowBest}. A method of its own, not inlined
		 * ({@code jvm.options}), so that the values {@link #rows} keeps across its rows are not live in the cell loop: on
		 * AArch64 C2 kept the loop's arrays and bestScore in stack slots there, reloaded every cell.
		 */
		private int row(GappedAligner g, int[] best, int[] bestGap, int[] mat, int mask, byte[] bs, byte[] tb, int base,
				int row, int ai, int first, int bSize, int bestScore, int xDrop) {
			final int goe = this.goe, ext = this.ext;
			int score = MININT, gapRow = MININT, last = first;
			for (int bi = first; bi < bSize; bi++) {
				int gapCol = bestGap[bi];
				int next = best[bi] + mat[(row + bs[bi]) & mask];
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
						g.aOffset = ai;
						g.bOffset = bi;
					}
					gapRow -= ext;
					gapCol -= ext;
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
			g.rowFirst = first;
			g.rowGapRow = gapRow;
			g.rowBest = bestScore;
			return last;
		}

		/**
		 * The rows of A after the initial row, whose band ends at {@code bSize}, then the traceback from the best cell
		 * to the start into {@code ops}; the best score.
		 */
		@Override
		public int rows(GappedAligner g, byte[] a, int aBase, int m, byte[] bs, int n, boolean reverse, int xDrop,
				int bSize, EditScript ops) {
			final int goe = this.goe, ext = this.ext;
			final int[] best = g.best, bestGap = g.bestGap, rowStart = g.rowStart;
			// The 32-wide matrix: with its index masked with mat.length - 1 (1024 entries), C2 drops the bounds check.
			final int[] mat = g.matrix.wide;
			final int mask = mat.length - 1;
			final int numExtraCells = ext > 0 ? xDrop / ext + 3 : n + 3;
			final int[] saved = g.savedBest, savedGap = g.savedGap;
			byte[] tb = g.trace;
			int used = bSize + 1;
			int bestScore = 0, first = 0;
			// Rows computed again with traceback, through this one, where the score-only rows found a better cell.
			int traceUntil = 0;
			for (int ai = 1; ai <= m; ai++) {
				// Room for this row: from 'first' to at most numExtraCells beyond the previous row's end.
				g.ensureTrace(used + 2 + (bSize - first) + numExtraCells + 1);
				tb = g.trace;
				final int base = used + 1 - first; // tb[base + bi] is cell (ai, bi)
				rowStart[ai] = base;
				final int origFirst = first;
				final int row = letterA(a, aBase, m, ai, reverse) << 5;
				final int last = row(g, best, bestGap, mat, mask, bs, tb, base, row, ai, first, bSize, bestScore, xDrop);
				first = g.rowFirst;
				bestScore = g.rowBest;
				int gapRow = g.rowGapRow;
				final int bi = Math.max(origFirst, bSize);
				if (first == bSize) {
					break;
				}
				if (last < bSize - 1) {
					bSize = last + 1;
				} else {
					while (gapRow >= bestScore - xDrop && bSize <= n) {
						best[bSize] = gapRow;
						bestGap[bSize] = gapRow - goe;
						gapRow -= ext;
						if (base + bSize >= tb.length) {
							g.ensureTrace(base + bSize + numExtraCells + 2);
							tb = g.trace;
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
				// Over half of the cells lie in rows after the best cell's, where only the scores matter: whether a better
				// cell comes. So AHEAD rows after the best, the next rows are computed score only, from a copy of the band.
				// Without a better cell (score > bestScore, as above) the traceback has its rows; with one, the band is
				// put back and the rows are computed again with traceback, up to that cell's row at least.
				if (ai >= g.aOffset + AHEAD && ai >= traceUntil && ai < m) {
					System.arraycopy(best, first, saved, first, bSize - first);
					System.arraycopy(bestGap, first, savedGap, first, bSize - first);
					final int better = scoreRows(g, a, aBase, m, bs, n, reverse, xDrop, ai + 1, first, bSize, bestScore,
							mat, mask);
					if (better < 0) {
						break;
					}
					System.arraycopy(saved, first, best, first, bSize - first);
					System.arraycopy(savedGap, first, bestGap, first, bSize - first);
					traceUntil = better;
				}
			}

			// Follow the traceback from the best cell to the start.
			int ai = g.aOffset, bi = g.bOffset;
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
	}

	/** The row loop of {@link #semiGapped}: {@link Rows}, or a copy of it with constant gap costs. */
	interface RowLoops {
		int rows(GappedAligner g, byte[] a, int aBase, int m, int n, boolean reverse, int xDrop, int bSize);
	}

	/**
	 * The row loops with the gap costs as fields: {@code goe} (gapOpen + gapExtend) and {@code ext} (gapExtend). The
	 * loops read them once into locals. {@link #of} gives a copy of this class whose reads of the two fields are the
	 * constants instead, made with the class-file API and defined as a hidden class in this nest (so it reads the
	 * aligner's private arrays as this class does). C2 then folds the gap costs into the loops for whatever costs the
	 * search uses.
	 */
	static final class Rows implements RowLoops {
		/** One copy per pair of gap costs. */
		private static final java.util.Map<Key, RowLoops> COPIES = new java.util.concurrent.ConcurrentHashMap<>();
		/** A native image cannot define classes at run time: there the fields stay. */
		static final boolean NATIVE_IMAGE = System.getProperty("org.graalvm.nativeimage.imagecode") != null;

		/** BLOSUM62 padded to 32 columns: row {@code a << 5}, column {@code b}. */
		static final int[] BLOSUM62 = new int[32 * 32];
		static {
			for (int a = 0; a < N; a++) {
				for (int b = 0; b < N; b++) {
					BLOSUM62[a << 5 | b] = Blosum62.score(a, b);
				}
			}
		}
		/**
		 * The score matrix padded to 32 columns, as a constant of a copy ({@link #of} with a matrix); null here, where
		 * the loop reads the aligner's matrix.
		 */
		private static final int[] MATRIX = null;

		private final int goe, ext;

		Rows(int goe, int ext) {
			this.goe = goe;
			this.ext = ext;
		}

		/**
		 * The row loop for these gap costs, and for {@code matrix} (padded as {@link #BLOSUM62}) unless it is null: a
		 * copy with them as constants, or this class where that fails.
		 */
		static RowLoops of(int goe, int ext, int[] matrix) {
			if (NATIVE_IMAGE) {
				return new Rows(goe, ext);
			}
			return COPIES.computeIfAbsent(new Key(goe, ext, matrix), key -> {
				try {
					return (RowLoops) copy(Rows.class, goe, ext, matrix);
				} catch (Throwable e) {
					return new Rows(goe, ext);
				}
			});
		}

		/** {@link GappedAligner#semiGapped}'s rows of A after the initial row, whose band ends at {@code bSize}. */
		public int rows(GappedAligner g, byte[] a, int aBase, int m, int n, boolean reverse, int xDrop, int bSize) {
			final int goe = this.goe, ext = this.ext;
			final int[] dp = g.dp;
			final int[] mat = g.matrix.scores, fixed = MATRIX;
			final byte[] bs = g.bLetters;
			// The scores of A's current letter. Allocated here with a constant length, so that C2 knows an index masked
			// with 31 is in bounds and does not check it every cell.
			final int[] rowScores = new int[32];
			int score;
			int bestScore = 0, first = 0;
			// bestScore - xDrop, kept up to date, so that the inner loop does not need xDrop.
			int floor = -xDrop;
			for (int ai = 1; ai <= m; ai++) {
				// With a constant matrix (a copy's), its row is read in place: with both indices masked C2 knows the
				// index is within its 1024 entries. Otherwise the row is copied, as above. C2 keeps one of the two.
				final int letter = letterA(a, aBase, m, ai, reverse);
				final int row = (letter & 31) << 5;
				if (fixed == null) {
					System.arraycopy(mat, letter * N, rowScores, 0, N);
				}
				score = MININT;
				int gapRow = MININT, last = first;
				for (int bi = first; bi < bSize; bi++) {
					int gapCol = dp[2 * bi + 1];
					final int b = bs[bi] & 31;
					int next = dp[2 * bi] + (fixed != null ? fixed[row | b] : rowScores[b]);
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
							g.aOffset = ai;
							g.bOffset = bi;
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
			return bestScore;
		}
	}
}
