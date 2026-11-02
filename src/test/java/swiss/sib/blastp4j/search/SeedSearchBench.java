package swiss.sib.blastp4j.search;

import java.nio.file.Path;
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

import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.seed.SeedParameters;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Step 5: the seeding stage (lookup table, two-hit trigger, ungapped extension) of one query against all of Swiss-Prot,
 * on one thread. The queries are HBA_HUMAN (142 residues, 6,339 subjects reach the gapped stage), CDK1_HUMAN (297,
 * 15,883) and BRCA1_HUMAN (1,863, 65,541). Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3)
@Measurement(iterations = 10)
@Fork(value = 1, jvmArgs = { "-Xmx2g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class SeedSearchBench {
	@Param({ "HBA_HUMAN", "CDK1_HUMAN", "BRCA1_HUMAN" })
	public String query;

	private SeedSearch search;
	private SeedParameters parameters;
	private byte[] sequence;
	private ExecutorService pool;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		SequenceDatabase db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		for (int i = 0; i < queries.size(); i++) {
			if (queries.id(i).endsWith("|" + query)) {
				sequence = queries.sequence(i);
			}
		}
		search = new SeedSearch(db);
		parameters = SeedParameters.blastp(sequence,
				new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues()), 10, BlastpSearch.MIN_SUBJECT_LENGTH);
		pool = Executors.newSingleThreadExecutor();
	}

	@TearDown
	public void tearDown() {
		pool.shutdown();
	}

	@Benchmark
	public List<SeedSearch.Seeds> seed() throws InterruptedException {
		return search.search(sequence, parameters, pool, 1);
	}
}
