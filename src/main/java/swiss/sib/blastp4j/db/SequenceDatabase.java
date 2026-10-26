package swiss.sib.blastp4j.db;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.List;

/**
 * Protein sequences laid out as in a BLAST {@code .psq} file: NCBIstdaa codes that start with a sentinel and have a
 * sentinel ({@link Alphabet#GAP}) after every sequence. The sentinels bound extensions without extra checks.
 * <p>
 * The residues are held in memory segments, one per volume: the memory-mapped {@code .psq} files of a BLAST database,
 * or a heap array for a FASTA file. A database is therefore not limited to the 2 GB of a Java array, and a mapped one
 * is paged in by the operating system rather than copied onto the Java heap. The search copies one subject at a time,
 * with its two sentinels, into a buffer of its own ({@link #withSentinels}).
 * <p>
 * Thread safe until {@link #close()}, which unmaps the files.
 */
public final class SequenceDatabase implements AutoCloseable {
	/** Offsets in a {@code .pin} file: big-endian 32-bit integers, not necessarily aligned. */
	static final ValueLayout.OfInt PIN_INT = ValueLayout.JAVA_INT_UNALIGNED.withOrder(ByteOrder.BIG_ENDIAN);

	/**
	 * One volume: its residues, the start of each of its sequences in them (unsigned 32-bit integers, one more than the
	 * number of sequences: the last is the end of the final sentinel), and the length of its longest sequence.
	 */
	record Volume(MemorySegment residues, MemorySegment starts, ValueLayout.OfInt layout, int size, int maxLength) {
		long start(int i) {
			return Integer.toUnsignedLong(starts.getAtIndex(layout, i));
		}
	}

	private final Arena arena;
	private final Volume[] volumes;
	/** Index of the first sequence of each volume, plus the total. */
	private final int[] firstOid;
	private final Headers headers;
	private final long totalResidues;
	private final int maxLength;

	/**
	 * A database in one heap array.
	 *
	 * @param residues the {@code .psq} layout: sentinel, sequence 0, sentinel, sequence 1, sentinel, ...
	 * @param starts   index in {@code residues} of the first residue of every sequence, plus one past the last sentinel
	 * @param headers  identifiers and titles
	 */
	SequenceDatabase(byte[] residues, int[] starts, Headers headers) {
		this(null, List.of(new Volume(MemorySegment.ofArray(residues), MemorySegment.ofArray(starts),
				ValueLayout.JAVA_INT, starts.length - 1, maxLength(starts))), headers);
	}

	private static int maxLength(int[] starts) {
		int max = 0;
		for (int i = 0; i + 1 < starts.length; i++) {
			max = Math.max(max, starts[i + 1] - starts[i] - 1);
		}
		return max;
	}

	/**
	 * @param arena   the arena the segments were mapped in, closed by {@link #close()}; null for heap segments
	 * @param volumes the volumes, in order
	 * @param headers identifiers and titles of all sequences
	 */
	SequenceDatabase(Arena arena, List<Volume> volumes, Headers headers) {
		this.arena = arena;
		this.volumes = volumes.toArray(Volume[]::new);
		this.firstOid = new int[this.volumes.length + 1];
		long total = 0;
		int max = 0;
		for (int v = 0; v < this.volumes.length; v++) {
			Volume vol = this.volumes[v];
			firstOid[v + 1] = Math.addExact(firstOid[v], vol.size);
			total += vol.start(vol.size) - 1 - vol.size;
			max = Math.max(max, vol.maxLength);
		}
		if (firstOid[this.volumes.length] != headers.size()) {
			throw new IllegalArgumentException("need one start per sequence plus the end");
		}
		this.headers = headers;
		this.totalResidues = total;
		this.maxLength = max;
	}

	/** Number of sequences. */
	public int size() {
		return headers.size();
	}

	/** Total number of residues, without sentinels (the database length of the statistics). */
	public long totalResidues() {
		return totalResidues;
	}

	/** Length of the longest sequence. */
	public int maxLength() {
		return maxLength;
	}

	/** The identifiers of sequence {@code i} as BLAST prints them as {@code sseqid}. */
	public String id(int i) {
		return headers.seqid(i);
	}

	/** The first word of the header of sequence {@code i} (for a query: its {@code qseqid}). */
	public String name(int i) {
		return headers.name(i);
	}

	/** {@code sacc} or {@code saccver} of sequence {@code i}. */
	public String accession(int i, boolean withVersion) {
		return headers.accession(i, withVersion);
	}

	/** {@code stitle}: the title of sequence {@code i}. */
	public String title(int i) {
		return headers.title(i);
	}

	/** The volume that holds sequence {@code i}. */
	private int volume(int i) {
		if (volumes.length == 1) {
			return 0;
		}
		int v = Arrays.binarySearch(firstOid, i);
		// An exact match can be the start of an empty volume; take the last volume that starts at i.
		if (v >= 0) {
			while (v + 1 < volumes.length && firstOid[v + 1] == i) {
				v++;
			}
			return v;
		}
		return -v - 2;
	}

	/** Length of sequence {@code i}. */
	public int length(int i) {
		int v = volume(i);
		Volume vol = volumes[v];
		int k = i - firstOid[v];
		return (int) (vol.start(k + 1) - vol.start(k) - 1);
	}

	/**
	 * Sequence {@code i} with the sentinels before and after it, as a view of the mapped residues: the residues are at
	 * {@code [1, 1 + length(i))}. Nothing is copied, and the codes are the database's (NCBIstdaa).
	 */
	public MemorySegment withSentinelsView(int i) {
		int v = volume(i);
		Volume vol = volumes[v];
		int k = i - firstOid[v];
		long start = vol.start(k);
		int n = (int) (vol.start(k + 1) - start) + 1;
		return vol.residues.asSlice(start - 1, n);
	}

	/**
	 * Copies sequence {@code i} with the sentinels before and after it into {@code buffer}: the residues are at
	 * {@code [1, 1 + length(i))}. Returns {@code buffer}, or a new array of {@code maxLength() + 2} bytes when it is
	 * null or too small.
	 */
	public byte[] withSentinels(int i, byte[] buffer) {
		int v = volume(i);
		Volume vol = volumes[v];
		int k = i - firstOid[v];
		long start = vol.start(k);
		int n = (int) (vol.start(k + 1) - start) + 1;
		if (buffer == null || buffer.length < n) {
			buffer = new byte[Math.max(n, maxLength + 2)];
		}
		MemorySegment.copy(vol.residues, ValueLayout.JAVA_BYTE, start - 1, buffer, 0, n);
		return buffer;
	}

	/** Copies sequence {@code i} into {@code dst} from {@code offset}, without sentinels; returns its length. */
	public int copy(int i, byte[] dst, int offset) {
		int v = volume(i);
		Volume vol = volumes[v];
		int k = i - firstOid[v];
		long start = vol.start(k);
		int n = (int) (vol.start(k + 1) - start - 1);
		MemorySegment.copy(vol.residues, ValueLayout.JAVA_BYTE, start, dst, offset, n);
		return n;
	}

	/** A copy of sequence {@code i} as NCBIstdaa codes. */
	public byte[] sequence(int i) {
		int v = volume(i);
		Volume vol = volumes[v];
		int k = i - firstOid[v];
		long start = vol.start(k);
		return vol.residues.asSlice(start, vol.start(k + 1) - start - 1).toArray(ValueLayout.JAVA_BYTE);
	}

	/** Sequence {@code i} as letters. */
	public String letters(int i) {
		StringBuilder sb = new StringBuilder(length(i));
		for (byte b : sequence(i)) {
			sb.append(Alphabet.decode(b));
		}
		return sb.toString();
	}

	/** Number of volumes. */
	public int volumes() {
		return volumes.length;
	}

	/** The residues of volume {@code v}, in the {@code .psq} layout; read only. */
	public MemorySegment residues(int v) {
		return volumes[v].residues.asReadOnly();
	}

	/** Unmaps the files of a BLAST database; the database must not be used afterwards. */
	@Override
	public void close() {
		if (arena != null) {
			arena.close();
		}
	}
}
