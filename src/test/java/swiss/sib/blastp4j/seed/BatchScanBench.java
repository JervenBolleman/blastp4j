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
 * The batch scan alone: {@link BatchWordFinder} over all of the database for a batch of queries, on one thread, as
 * blastp's preliminary stage scans one batch. This is the shape the full search spends the scan in, which
 * {@link WordFinderBench} (one query, {@link WordFinder}) is not: the batch's lookup table is the one whose cells and
 * entries the scan's inner loops read, and its size is what decides how many hits a subject position has.
 * <p>
 * Why it exists: the scan is where the search spends its seeding, and a full search dilutes it. Halving blastp's batch
 * size, which doubles the number of scans, costs 2.95% of a full search's instructions, so the fixed cost per subject
 * position is about 160 instructions of the roughly 1660 a position costs at a batch of 10000 residues -- the rest
 * scales with the hits, which the batch's size decides.
 * <p>
 * What it cannot do: resolve a change of a few instructions per position. On an i7-6700HQ its instruction count
 * repeats to about 1.5% (standard deviation over three runs at {@code batchSize=1000}, unforked, {@code -wi 2 -i 5}),
 * and its {@code ms/op} to about 3%, both worse than a full search's 0.9% on instructions, because only about 28 G
 * instructions run and the JIT's and the allocator's variation are a large share of that. For a change worth a few
 * instructions per position, count them in the compiled loop body instead ({@code -XX:CompileCommand=print,
 * swiss.sib.blastp4j.seed.BatchWordFinder*::find}; the class is hidden, so the pattern needs the wildcard) and use
 * this benchmark for changes large enough to clear its spread.
 * <p>
 * The output is the number of ungapped alignments found, which must not change when the code is made faster. Needs the
 * data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 1, jvmArgs = { "-Xmx2g", "--add-modules", "jdk.incubator.vector" })
@State(Scope.Benchmark)
public class BatchScanBench {
	/**
	 * The batch's size in residues. 10000 is blastp's own default ({@code BLASTP_BATCH_SIZE}), the shape the search
	 * runs in; 1000 is there to make the scan's fixed cost per subject position visible. At 10000 a position costs
	 * about 1660 instructions, only about 160 of them fixed, so a change of a few instructions per position is lost
	 * in the hit work; a batch a tenth of the size has about a tenth of the hits per position, which puts the fixed
	 * cost at roughly half and makes such a change a few percent of what is measured.
	 */
	@Param({ "1000", "10000" })
	public int batchSize;

	private SequenceDatabase db;
	private byte[][] queries;
	private LookupTable[] tables;
	private SeedParameters[] seeds;
	private BatchWordFinder.Table table;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		db = BlastDbReader.exists(dir.resolve("blastdb/sprot")) ? BlastDbReader.read(dir.resolve("blastdb/sprot"))
				: FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase qs = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		BlastpSearch search = new BlastpSearch(db, statistics, 10, false);
		// As BlastpSearch.batch does: queries until the batch holds batchSize residues, skipping those whose
		// composition gives no Karlin-Altschul parameters.
		List<byte[]> batch = new ArrayList<>();
		List<SeedParameters> parameters = new ArrayList<>();
		long residues = 0;
		for (int i = 0; i < qs.size() && residues < batchSize; i++) {
			try {
				parameters.add(search.seedParameters(qs.sequence(i)));
				batch.add(qs.sequence(i));
				residues += qs.length(i);
			} catch (IllegalArgumentException e) {
				// Skipped, as blastp skips it.
			}
		}
		if (batch.size() < 2) {
			throw new IllegalStateException("a batch of " + batch.size() + " queries: the scan needs two or more");
		}
		queries = batch.toArray(byte[][]::new);
		seeds = parameters.toArray(SeedParameters[]::new);
		tables = new LookupTable[queries.length];
		for (int i = 0; i < queries.length; i++) {
			tables[i] = LookupTable.of(queries[i], seeds[i].threshold());
		}
		table = BatchWordFinder.Table.of(queries, tables, seeds[0].window());
	}

	@TearDown
	public void tearDown() {
		db.close();
	}

	@Benchmark
	@SuppressWarnings("unchecked")
	public long scan(Blackhole bh) {
		// A fresh scan, as each batch of the search gets: BatchWordFinder.of generates a copy of the class with the
		// batch's constants, which is once per batch against a whole database scan.
		BatchScan finder = BatchWordFinder.of(queries, tables, seeds, table);
		final List<UngappedHsp>[] out = new List[queries.length];
		for (int i = 0; i < out.length; i++) {
			out[i] = new ArrayList<>();
		}
		long count = 0;
		byte[] subject = null;
		for (int i = 0; i < db.size(); i++) {
			subject = db.withSentinels(i, subject);
			finder.find(subject, 1, db.length(i), out);
			for (int q = 0; q < out.length; q++) {
				count += out[q].size();
				out[q].clear();
			}
		}
		bh.consume(count);
		return count;
	}
}
