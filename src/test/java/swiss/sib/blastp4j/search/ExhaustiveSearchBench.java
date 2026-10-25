package swiss.sib.blastp4j.search;

import java.nio.file.Path;
import java.util.List;
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
import org.openjdk.jmh.annotations.Warmup;

import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Step 4: one query against all of Swiss-Prot (209 million residues) with the exhaustive Smith-Waterman search, on
 * one thread. Single-shot time per search; GCUPS = query length x 209,017,843 / time. The queries are HBA_HUMAN (142
 * residues) and CDK1_HUMAN (297). Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.SECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 1, jvmArgs = { "-Xmx2g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class ExhaustiveSearchBench {
	@Param({ "HBA_HUMAN", "CDK1_HUMAN" })
	public String query;

	private ExhaustiveSearch search;
	private byte[] sequence;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		SequenceDatabase db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		search = new ExhaustiveSearch(db, new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues()), 10);
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		for (int i = 0; i < queries.size(); i++) {
			if (queries.id(i).endsWith("|" + query)) {
				sequence = queries.sequence(i);
			}
		}
	}

	@Benchmark
	public List<ExhaustiveSearch.Hit> search() {
		return search.search(sequence);
	}
}
