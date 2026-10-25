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

import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Steps 6 and 7: the complete search (seeding, preliminary gapped extensions, traceback, hit list) of one query
 * against all of Swiss-Prot, on one thread, without composition-based statistics (cbs0) and with blastp's default
 * mode 2 (cbs2). The queries are HBA_HUMAN (142 residues),
 * CDK1_HUMAN (297) and BRCA1_HUMAN (1,863). Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3)
@Measurement(iterations = 10)
@Fork(value = 1, jvmArgs = "-Xmx2g")
@State(Scope.Benchmark)
public class BlastpSearchBench {
	@Param({ "HBA_HUMAN", "CDK1_HUMAN", "BRCA1_HUMAN" })
	public String query;

	/** cbs0: without composition-based statistics; cbs2: blastp's default mode 2. */
	@Param({ "cbs0", "cbs2" })
	public String mode;

	private BlastpSearch search;
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
		search = new BlastpSearch(db, new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues()), 10,
				mode.equals("cbs2"));
		pool = Executors.newSingleThreadExecutor();
	}

	@TearDown
	public void tearDown() {
		pool.shutdown();
	}

	@Benchmark
	public List<BlastpSearch.Hit> search() throws InterruptedException {
		return search.search(sequence, pool, 1);
	}
}
