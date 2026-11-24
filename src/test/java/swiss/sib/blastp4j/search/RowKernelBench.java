package swiss.sib.blastp4j.search;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import swiss.sib.blastp4j.align.GappedAligner;
import swiss.sib.blastp4j.db.BlastDbReader;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.seed.UngappedHsp;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * The score-only X-drop extension alone ({@code Blast_SemiGappedAlign}: {@link GappedAligner#extendScoreOnly}, so its
 * prologue, the copy of B's letters and the row loop), replayed over a corpus of the real extensions of one query
 * against all of Swiss-Prot. Everything around the kernel — the scan, the two-hit logic, the sort of the seeds, the
 * containment tree, the HSP objects — is done once in the setup and is not measured, so that a change to the row loop
 * shows here at its full size rather than as a few percent of a search.
 * <p>
 * This is step 8's benchmark of the row kernel alone (see {@code docs/acceleration-plan.md}): on AArch64 the lanes of
 * {@code SemiGappedLanes} are off, so these extensions run scalar, and a kernel that vectorises within one row has to
 * be measured against this before it is worth wiring into the search.
 * <p>
 * The corpus is the best seed of every subject that seeding hit, which {@code BlastpSearch.gappedScores} always
 * extends: its sort puts that seed first and the containment tree is still empty there. The later seeds of a subject
 * are left out because the tree skips most of them, so including them would weight the corpus by something this
 * benchmark does not model. The start point is the real one ({@link GappedAligner#startForGapped}) and the X-drop is
 * BLAST's for these scoring parameters.
 * <p>
 * The output is the sum of the scores, which <b>must not change</b> when the kernel is made faster: it is the
 * correctness check that goes with the timing. Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3)
@Measurement(iterations = 10)
@Fork(value = 1, jvmArgs = { "-Xmx3g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class RowKernelBench {
	/** BLAST's X-drop for the preliminary gapped stage is 15 bits, as {@code BlastpSearch}. */
	private static final double LN2 = 0.69314718055994530941723212145818;

	@Param({ "CDK1_HUMAN", "BRCA1_HUMAN" })
	public String query;

	/** Whether to measure the restricted extension ({@code s_RestrictedGappedAlign}) instead of the full one. */
	@Param({ "false" })
	public boolean restricted;

	/** One extension of the corpus: the subject as the stage sees it, its length, and the start point. */
	public record Extension(byte[] residues, int length, int qOff, int sOff) {
	}

	private byte[] sequence;
	private List<Extension> corpus;
	private GappedAligner aligner;
	private int xDrop;
	/** The sum of the scores of the last invocation, printed at the end as the correctness check. */
	private long scores;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		try (SequenceDatabase db = BlastDbReader.exists(dir.resolve("blastdb/sprot"))
				? BlastDbReader.read(dir.resolve("blastdb/sprot"))
				: FastaReader.read(dir.resolve("uniprot_sprot.fasta"))) {
			SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
			for (int i = 0; i < queries.size(); i++) {
				if (queries.id(i).endsWith("|" + query)) {
					sequence = queries.sequence(i);
				}
			}
			if (sequence == null) {
				throw new IllegalArgumentException("no query " + query);
			}
			Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
			// As BlastpSearch computes it from the scoring parameters' lambda, which is all this needs of the search.
			xDrop = (int) (15 * LN2 / statistics.parameters().lambda());
			BlastpSearch search = new BlastpSearch(db, statistics, 10, false);
			ExecutorService pool = Executors.newSingleThreadExecutor();
			List<SeedSearch.Seeds> seeds;
			try {
				seeds = new SeedSearch(db).search(sequence, search.seedParameters(sequence), pool, 1);
			} finally {
				pool.shutdown();
			}
			aligner = new GappedAligner(11, 1);
			corpus = new ArrayList<>(seeds.size());
			long residues = 0, cells = 0;
			for (SeedSearch.Seeds s : seeds) {
				// A copy per subject, not the one reused buffer of the search: the corpus outlives the loop.
				final byte[] subject = db.withSentinels(s.subject(), null);
				final int length = db.length(s.subject());
				List<UngappedHsp> hsps = new ArrayList<>(s.hsps());
				hsps.sort(BlastpSearch.SEED_ORDER);
				final UngappedHsp seed = hsps.get(0);
				final int qOff = aligner.startForGapped(sequence, subject, 1, seed.queryStart(), seed.length(),
						seed.subjectStart(), seed.length());
				final int sOff = seed.subjectSeed() + qOff - seed.querySeed();
				corpus.add(new Extension(subject, length, qOff, sOff));
				residues += length;
				cells += (long) sequence.length * length;
			}
			System.out.printf("%ncorpus: %d extensions, %d subject residues, query %d residues, %.1f G cells if full%n",
					corpus.size(), residues, sequence.length, cells / 1e9);
		}
	}

	@TearDown
	public void tearDown() {
		// Printed rather than asserted: the sum belongs to the query, and what matters is that it does not move
		// between two kernels measured on the same corpus.
		System.out.printf("%nscores: %d over %d extensions of %s%n", scores, corpus.size(), query);
		// The subjects are the corpus's own copies, so the database is already closed.
		corpus = null;
	}

	/**
	 * The extensions of the corpus, in order, through one aligner, as the gapped stage runs them on one thread. The
	 * aligner is the setup's, so its {@code dp} has grown already and the timing is the kernel's, not the first row's.
	 */
	@Benchmark
	public long extend(Blackhole bh) {
		long sum = 0;
		for (int i = 0, n = corpus.size(); i < n; i++) {
			final Extension e = corpus.get(i);
			aligner.extendScoreOnly(sequence, e.residues(), 1, e.length(), e.qOff(), e.sOff(), restricted, xDrop);
			sum += aligner.score;
		}
		bh.consume(sum);
		scores = sum;
		return sum;
	}
}
