package swiss.sib.blastp4j.compo;

import java.nio.file.Path;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

import swiss.sib.blastp4j.db.BlastDbReader;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;

/**
 * SEG with the composition-based statistics parameters over all of Swiss-Prot, one thread: the segments of every
 * subject. The output is the number of segments and the number of masked letters, which must not change when the code
 * is made faster. Needs the data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 1, jvmArgs = "-Xmx2g")
@State(Scope.Benchmark)
public class SegBench {
	private SequenceDatabase db;
	private Seg seg;

	@Setup
	public void setup() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		db = BlastDbReader.exists(dir.resolve("blastdb/sprot")) ? BlastDbReader.read(dir.resolve("blastdb/sprot"))
				: FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		seg = Seg.forCompositionBasedStatistics();
	}

	@TearDown
	public void tearDown() {
		db.close();
	}

	@Benchmark
	public long segments(Blackhole bh) {
		long count = 0;
		byte[] subject = null;
		for (int i = 0; i < db.size(); i++) {
			subject = db.withSentinels(i, subject);
			for (int[] s : seg.segments(subject, 1, db.length(i))) {
				count += 1 + s[1] - s[0];
			}
		}
		bh.consume(count);
		return count;
	}
}
