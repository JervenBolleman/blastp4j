package swiss.sib.blastp4j.search;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import swiss.sib.blastp4j.db.Alphabet;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.stats.Blosum62;
import swiss.sib.blastp4j.stats.Statistics;
import swiss.sib.cssw4j.Aligner;
import swiss.sib.cssw4j.Alignment;
import swiss.sib.cssw4j.Profile;
import swiss.sib.cssw4j.ScoringMatrix;

/**
 * The exhaustive baseline: the optimal local alignment (Smith-Waterman, cssw4j) of the query against every subject of
 * the database, with BLAST's scoring system (its BLOSUM62 and gap costs 11/1) and its statistics. Reports one hit per
 * subject, the best alignment, when its E-value is at most the cutoff.
 * <p>
 * BLAST charges {@code open + k * extend} for a gap of k residues, SSW {@code open + (k - 1) * extend}, so BLAST's
 * 11/1 is SSW's 12/1. Without composition-based statistics, every BLAST HSP is a local alignment under the same
 * scores, so the score found here is an upper bound for BLAST's best HSP with each subject.
 * <p>
 * Instances are immutable; {@link #search} may be called from several threads.
 */
public final class ExhaustiveSearch {
	/** One subject with its best local alignment; positions are 0-based and inclusive. */
	public record Hit(int subject, int score, double evalue, int queryStart, int queryEnd, int subjectStart,
			int subjectEnd) {
	}

	/** BLAST's BLOSUM62 for cssw4j: NCBIstdaa codes are the matrix indices, rows are subject letters. */
	public static final ScoringMatrix MATRIX;

	static {
		byte[] m = new byte[Blosum62.SIZE * Blosum62.SIZE];
		for (int s = 0; s < Blosum62.SIZE; s++) {
			for (int q = 0; q < Blosum62.SIZE; q++) {
				m[s * Blosum62.SIZE + q] = (byte) Blosum62.score(q, s);
			}
		}
		byte[] table = new byte[128];
		for (int c = 0; c < 128; c++) {
			table[c] = Alphabet.encode(c);
		}
		MATRIX = ScoringMatrix.of(Blosum62.SIZE, m, table);
	}

	/** Subjects per work unit. */
	private static final int CHUNK = 2048;

	private final SequenceDatabase db;
	private final Statistics statistics;
	private final double evalueCutoff;
	private final int gapOpen, gapExtend;
	private final int maxLength;

	/**
	 * @param db           the subjects
	 * @param statistics   E-values; its scoring parameters give the gap costs (BLAST convention)
	 * @param evalueCutoff report subjects with an E-value at most this
	 */
	public ExhaustiveSearch(SequenceDatabase db, Statistics statistics, double evalueCutoff) {
		this.db = db;
		this.statistics = statistics;
		this.evalueCutoff = evalueCutoff;
		this.gapOpen = statistics.parameters().gapOpen() + statistics.parameters().gapExtend();
		this.gapExtend = statistics.parameters().gapExtend();
		int max = 0;
		for (int i = 0; i < db.size(); i++) {
			max = Math.max(max, db.length(i));
		}
		this.maxLength = max;
	}

	/**
	 * Aligns a query (NCBIstdaa codes) against all subjects, using {@code pool} with {@code tasks} tasks. Hits are
	 * ordered by E-value, then score (descending), then subject index.
	 */
	public List<Hit> search(byte[] query, ExecutorService pool, int tasks) throws InterruptedException {
		Profile profile = Profile.of(query, query.length, MATRIX, Profile.ScoreSize.BOTH);
		AtomicInteger next = new AtomicInteger();
		List<Future<List<Hit>>> parts = new ArrayList<>();
		for (int t = 0; t < tasks; t++) {
			parts.add(pool.submit(() -> scan(profile, next)));
		}
		List<Hit> hits = new ArrayList<>();
		try {
			for (Future<List<Hit>> part : parts) {
				hits.addAll(part.get());
			}
		} catch (ExecutionException e) {
			throw new IllegalStateException(e.getCause());
		}
		hits.sort(Comparator.comparingDouble(Hit::evalue).thenComparing(Comparator.comparingInt(Hit::score).reversed())
				.thenComparingInt(Hit::subject));
		return hits;
	}

	/** Single-threaded {@link #search}. */
	public List<Hit> search(byte[] query) {
		try (ExecutorService pool = Executors.newSingleThreadExecutor()) {
			return search(query, pool, 1);
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IllegalStateException(e);
		}
	}

	/** Takes chunks of subjects from {@code next} until none are left. */
	private List<Hit> scan(Profile profile, AtomicInteger next) {
		final Aligner aligner = new Aligner();
		final byte[] residues = db.residues();
		final byte[] subject = new byte[maxLength];
		final int queryLength = profile.length();
		final List<Hit> hits = new ArrayList<>();
		for (int chunk; (chunk = next.getAndIncrement()) * CHUNK < db.size();) {
			int end = Math.min(db.size(), (chunk + 1) * CHUNK);
			for (int i = chunk * CHUNK; i < end; i++) {
				int length = db.length(i);
				System.arraycopy(residues, db.start(i), subject, 0, length);
				Alignment a = aligner.align(profile, subject, length, gapOpen, gapExtend, 0, 0, 0, 0);
				if (a.score1() <= 0) {
					continue;
				}
				double evalue = statistics.evalue(a.score1(), queryLength, length);
				if (evalue <= evalueCutoff) {
					// Begin positions only for the few hits: one more (reverse) alignment.
					Alignment b = aligner.align(profile, subject, length, gapOpen, gapExtend, Aligner.FLAG_BEGIN, 0, 0,
							0);
					hits.add(new Hit(i, a.score1(), evalue, b.readBegin1(), b.readEnd1(), b.refBegin1(), b.refEnd1()));
				}
			}
		}
		return hits;
	}
}
