package swiss.sib.blastp4j.stats;

import java.util.List;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OperationsPerInvocation;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import swiss.sib.blastp4j.ReferenceData;
import swiss.sib.blastp4j.output.TabularFormat;

/**
 * Step 2: cost of the statistics for one reported alignment. Spouge's finite-size corrected E-value (two erfc calls)
 * against the classic {@code K m n e^(-lambda S)}, and the tabular formatting, over the 15,107 hits of the reference.
 * Average time per hit.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(2)
@State(Scope.Benchmark)
public class StatisticsBench {
	static final int HITS = 15_107;

	private int[] scores, queryLengths, subjectLengths;
	private Statistics stats;

	@Setup(Level.Trial)
	public void setup() throws Exception {
		List<String[]> rows = ReferenceData.load("human_vs_sprot.cbs0.tsv.gz");
		scores = new int[HITS];
		queryLengths = new int[HITS];
		subjectLengths = new int[HITS];
		for (int i = 0; i < HITS; i++) {
			String[] r = rows.get(i);
			scores[i] = Integer.parseInt(r[ReferenceData.SCORE]);
			queryLengths[i] = Integer.parseInt(r[ReferenceData.QLEN]);
			subjectLengths[i] = Integer.parseInt(r[ReferenceData.SLEN]);
		}
		stats = new Statistics(ScoringParameters.BLOSUM62_11_1, ReferenceData.SPROT_RESIDUES);
	}

	@Benchmark
	@OperationsPerInvocation(HITS)
	public void spougeEvalue(Blackhole bh) {
		for (int i = 0; i < HITS; i++) {
			bh.consume(stats.evalue(scores[i], queryLengths[i], subjectLengths[i]));
		}
	}

	@Benchmark
	@OperationsPerInvocation(HITS)
	public void classicEvalue(Blackhole bh) {
		ScoringParameters p = stats.parameters();
		double n = ReferenceData.SPROT_RESIDUES;
		for (int i = 0; i < HITS; i++) {
			bh.consume(p.k() * queryLengths[i] * n * Math.exp(-p.lambda() * scores[i]));
		}
	}

	@Benchmark
	@OperationsPerInvocation(HITS)
	public void spougeAndFormat(Blackhole bh) {
		for (int i = 0; i < HITS; i++) {
			bh.consume(TabularFormat.evalue(stats.evalue(scores[i], queryLengths[i], subjectLengths[i])));
			bh.consume(TabularFormat.bitScore(stats.bitScore(scores[i])));
		}
	}
}
