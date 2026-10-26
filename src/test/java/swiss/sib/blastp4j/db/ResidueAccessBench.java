package swiss.sib.blastp4j.db;

import java.io.IOException;
import java.io.OutputStream;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;

/**
 * The cost of reading residues through a memory segment against a byte array, in the word finder's inner loop. A 5-bit
 * rolling word over {@value #SIZE} residues of a file mapped in a shared arena; for each word, its neighbours from a
 * table, each with a branch and a write to a diagonal array when it passes. The same loop runs over the mapped segment
 * ({@code segment}), over a heap copy of the file ({@code array}), and over chunks copied from the segment into a small
 * heap buffer ({@code chunks}). Measured: segment about 16% slower than array; chunks as fast as array. See the commit
 * that reverted the zero-copy seeding.
 */
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(value = 1, jvmArgs = "-Xmx1g")
@State(Scope.Benchmark)
public class ResidueAccessBench {
	/** Residues per scan: about the size of Swiss-Prot. */
	static final int SIZE = 1 << 26;
	private static final int CHUNK = 4096;

	private Path file;
	private Arena arena;
	private MemorySegment segment;
	private byte[] array;
	private final byte[] buffer = new byte[CHUNK];
	/** Neighbours of each word, as the lookup table has them: {@code offsets[start[w] .. start[w + 1])}. */
	private int[] start, offsets;
	/** The last hit on each diagonal, written for the neighbours that score high enough. */
	private final int[] lastHit = new int[1 << 15];

	@Setup(Level.Trial)
	public void setup() throws IOException {
		file = Files.createTempFile("residues", ".bin");
		byte[] residues = new byte[SIZE];
		new Random(1).nextBytes(residues);
		for (int i = 0; i < SIZE; i++) {
			residues[i] = (byte) ((residues[i] & 0xff) % 28);
		}
		try (OutputStream out = Files.newOutputStream(file)) {
			out.write(residues);
		}
		array = residues;
		Random random = new Random(2);
		start = new int[(1 << 15) + 1];
		for (int w = 0; w < 1 << 15; w++) {
			start[w + 1] = start[w] + random.nextInt(3);
		}
		offsets = new int[start[1 << 15]];
		for (int k = 0; k < offsets.length; k++) {
			offsets[k] = random.nextInt();
		}
		arena = Arena.ofShared();
		try (FileChannel ch = FileChannel.open(file, StandardOpenOption.READ)) {
			segment = ch.map(FileChannel.MapMode.READ_ONLY, 0, ch.size(), arena);
		}
	}

	@TearDown(Level.Trial)
	public void tearDown() throws IOException {
		arena.close();
		Files.deleteIfExists(file);
	}

	/** One neighbour: a branch on its value, and a write to the diagonal's last hit when it passes. */
	private int hit(int neighbour, int position) {
		int diagonal = neighbour & 0x7fff;
		if ((neighbour & 3) != 0) {
			lastHit[diagonal] = position;
			return 1;
		}
		return 0;
	}

	@Benchmark
	public long array(Blackhole bh) {
		long count = 0;
		int w = 0;
		for (int i = 0; i < SIZE; i++) {
			w = ((w << 5) | array[i]) & 0x7fff;
			for (int k = start[w], end = start[w + 1]; k < end; k++) {
				count += hit(offsets[k], (int) i);
			}
		}
		bh.consume(count);
		return count;
	}

	@Benchmark
	public long segment(Blackhole bh) {
		long count = 0;
		int w = 0;
		for (long i = 0; i < SIZE; i++) {
			w = ((w << 5) | segment.get(ValueLayout.JAVA_BYTE, i)) & 0x7fff;
			for (int k = start[w], end = start[w + 1]; k < end; k++) {
				count += hit(offsets[k], (int) i);
			}
		}
		bh.consume(count);
		return count;
	}

	@Benchmark
	public long chunks(Blackhole bh) {
		long count = 0;
		int w = 0;
		for (int from = 0; from < SIZE; from += CHUNK) {
			MemorySegment.copy(segment, ValueLayout.JAVA_BYTE, from, buffer, 0, CHUNK);
			for (int i = 0; i < CHUNK; i++) {
				w = ((w << 5) | buffer[i]) & 0x7fff;
				for (int k = start[w], end = start[w + 1]; k < end; k++) {
					count += hit(offsets[k], i);
				}
			}
		}
		bh.consume(count);
		return count;
	}
}
