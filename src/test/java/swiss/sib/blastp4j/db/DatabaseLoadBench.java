package swiss.sib.blastp4j.db;

import java.nio.file.Path;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/**
 * Step 3: time to read UniProtKB/Swiss-Prot (575,748 sequences, 209 million residues) into memory, from the plain and
 * the gzip-compressed FASTA, and from the BLAST database made by makeblastdb. Single-shot time per load; needs the
 * data in {@code ~/data/uniprot} (or {@code -Dblastp4j.sprot}).
 */
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 2)
@Measurement(iterations = 5)
@Fork(value = 2, jvmArgs = "-Xmx2g")
@State(Scope.Benchmark)
public class DatabaseLoadBench {
	@Param({ "uniprot_sprot.fasta", "uniprot_sprot.fasta.gz", "blastdb/sprot" })
	public String file;

	@Benchmark
	public SequenceDatabase load() throws Exception {
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		Path path = dir.resolve(file);
		return BlastDbReader.exists(path) ? BlastDbReader.read(path) : FastaReader.read(path);
	}
}
