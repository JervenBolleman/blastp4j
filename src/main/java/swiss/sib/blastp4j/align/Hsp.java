package swiss.sib.blastp4j.align;

import java.util.Comparator;

/**
 * A gapped alignment of the query with one subject ({@code BlastHSP} of BLAST). Offsets are 0-based; the ends are
 * exclusive, as in BLAST.
 */
public final class Hsp {
	public int score;
	public int queryStart, queryEnd, subjectStart, subjectEnd;
	/** The start point of the gapped extension; the traceback starts from it again. */
	public int queryGappedStart, subjectGappedStart;
	public double evalue;
	/** The alignment operations, after the traceback; null before. */
	public EditScript editScript;

	public Hsp(int score, int queryStart, int queryEnd, int subjectStart, int subjectEnd, int queryGappedStart,
			int subjectGappedStart) {
		this.score = score;
		this.queryStart = queryStart;
		this.queryEnd = queryEnd;
		this.subjectStart = subjectStart;
		this.subjectEnd = subjectEnd;
		this.queryGappedStart = queryGappedStart;
		this.subjectGappedStart = subjectGappedStart;
	}

	/**
	 * BLAST's {@code ScoreCompareHSPs}: higher score first, then lower subject start, higher subject end, lower query
	 * start, higher query end.
	 */
	public static final Comparator<Hsp> BY_SCORE = (h1, h2) -> {
		int r = Integer.compare(h2.score, h1.score);
		if (r == 0) {
			r = Integer.compare(h1.subjectStart, h2.subjectStart);
		}
		if (r == 0) {
			r = Integer.compare(h2.subjectEnd, h1.subjectEnd);
		}
		if (r == 0) {
			r = Integer.compare(h1.queryStart, h2.queryStart);
		}
		if (r == 0) {
			r = Integer.compare(h2.queryEnd, h1.queryEnd);
		}
		return r;
	};

	/**
	 * {@code s_HSPIsContained} (with {@code CONTAINED_IN_HSP}): whether this alignment, scoring no more than
	 * {@code other}, has both its start and its end inside the query and subject ranges of {@code other}.
	 */
	public boolean isContainedIn(int score, Hsp other) {
		return isContained(score, queryStart, queryEnd, subjectStart, subjectEnd, other);
	}

	static boolean isContained(int score, int qStart, int qEnd, int sStart, int sEnd, Hsp other) {
		return score <= other.score && inside(other.queryStart, other.queryEnd, qStart, other.subjectStart,
				other.subjectEnd, sStart) && inside(other.queryStart, other.queryEnd, qEnd, other.subjectStart,
						other.subjectEnd, sEnd);
	}

	private static boolean inside(int qFrom, int qTo, int q, int sFrom, int sTo, int s) {
		return qFrom <= q && qTo >= q && sFrom <= s && sTo >= s;
	}

	@Override
	public String toString() {
		return "Hsp[score=" + score + ", q=" + queryStart + ".." + queryEnd + ", s=" + subjectStart + ".." + subjectEnd
				+ ", start=" + queryGappedStart + "/" + subjectGappedStart + ", " + editScript + "]";
	}
}
