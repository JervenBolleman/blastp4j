package swiss.sib.blastp4j.seed;

import java.nio.file.Path;
import java.util.ArrayList;
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
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import swiss.sib.blastp4j.db.BlastDbReader;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.search.BlastpSearch;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * The seeding stage alone: {@link WordFinder} over all of Swiss-Prot for one query, on one thread, as blastp's word
 * finding does it (word size 3, threshold 11, window 40, the query's ungapped X-drop). The output is the number of
 * ungapped alignments found, which must not change when the code is made faster. Needs the data in {@code ~/data/uniprot}
 * (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 1, jvmArgs = { "-Xmx2g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class WordFinderBench {
	@Param({ "HBA_HUMAN", "CDK1_HUMAN" })
	public String query;

	private SequenceDatabase db;
	private byte[] sequence;
	private LookupTable table;
	private SeedParameters seeds;

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
		seeds = new BlastpSearch(db, statistics, 10, false).seedParameters(sequence);
		table = LookupTable.of(sequence, seeds.threshold());
	}

	@TearDown
	public void tearDown() {
		db.close();
	}

	@Benchmark
	public long find(Blackhole bh) {
		WordScan finder = WordFinder.of(sequence, table, seeds);
		List<UngappedHsp> out = new ArrayList<>();
		long count = 0;
		byte[] subject = null;
		for (int i = 0; i < db.size(); i++) {
			subject = db.withSentinels(i, subject);
			if (finder.vectorScan()) {
				finder.findBlocks(subject, 1, db.length(i), out);
			} else {
				finder.find(subject, 1, db.length(i), out);
			}
			count += out.size();
			out.clear();
		}
		bh.consume(count);
		return count;
	}
}
