package swiss.sib.blastp4j.search;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
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
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.BlastDbReader;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.seed.UngappedHsp;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * The preliminary gapped stage alone ({@code BlastpSearch.gappedScores}: the score-only X-drop extensions, the
 * restricted alignment where it applies) over the seeds of all of Swiss-Prot for one query, on one thread. The seeds
 * are found once in the setup. The output is the number of gapped alignments, which must not change when the code is
 * made faster. Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 1, jvmArgs = { "-Xmx2g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class GappedScoreBench {
	@Param({ "CDK1_HUMAN", "BRCA1_HUMAN" })
	public String query;

	private SequenceDatabase db;
	private BlastpSearch search;
	private byte[] sequence;
	private List<SeedSearch.Seeds> seeds;
	private int cutoff;
	private ExecutorService pool;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		db = BlastDbReader.exists(dir.resolve("blastdb/sprot")) ? BlastDbReader.read(dir.resolve("blastdb/sprot"))
				: FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		for (int i = 0; i < queries.size(); i++) {
			if (queries.id(i).endsWith("|" + query)) {
				sequence = queries.sequence(i);
			}
		}
		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		search = new BlastpSearch(db, statistics, 10, false);
		cutoff = search.cutoff(sequence);
		pool = Executors.newSingleThreadExecutor();
		seeds = new SeedSearch(db).search(sequence, search.seedParameters(sequence), pool, 1);
	}

	@TearDown
	public void tearDown() {
		pool.shutdown();
		db.close();
	}

	@Benchmark
	public long gapped(Blackhole bh) {
		GappedAligner aligner = new GappedAligner(11, 1);
		long count = 0;
		byte[] subject = null;
		for (SeedSearch.Seeds s : seeds) {
			subject = db.withSentinels(s.subject(), subject);
			List<UngappedHsp> hsps = new ArrayList<>(s.hsps());
			List<Hsp> out = search.gappedScores(sequence, aligner, subject, db.length(s.subject()), hsps, cutoff);
			count += out.size();
		}
		bh.consume(count);
		return count;
	}
}
