package swiss.sib.blastp4j.db;

import java.util.List;

/**
 * Protein sequences in memory, laid out as in a BLAST {@code .psq} file: one array of NCBIstdaa codes that starts with
 * a sentinel and has a sentinel ({@link Alphabet#GAP}) after every sequence. The sentinels bound extensions without
 * extra checks.
 * <p>
 * Immutable and thread safe. The residue array is shared, not copied: callers must not modify it.
 */
public final class SequenceDatabase {
	private final byte[] residues;
	private final int[] starts;
	private final List<String> ids;
	private final long totalResidues;

	/**
	 * @param residues the {@code .psq} layout: sentinel, sequence 0, sentinel, sequence 1, sentinel, ...
	 * @param starts   index in {@code residues} of the first residue of every sequence, plus one past the last sentinel
	 * @param ids      sequence ids as BLAST prints them (the first word of the FASTA header)
	 */
	SequenceDatabase(byte[] residues, int[] starts, List<String> ids) {
		if (starts.length != ids.size() + 1) {
			throw new IllegalArgumentException("need one start per sequence plus the end");
		}
		this.residues = residues;
		this.starts = starts;
		this.ids = List.copyOf(ids);
		this.totalResidues = (long) starts[starts.length - 1] - 1 - ids.size();
	}

	/** Number of sequences. */
	public int size() {
		return ids.size();
	}

	/** Total number of residues, without sentinels (the database length of the statistics). */
	public long totalResidues() {
		return totalResidues;
	}

	/** The id of sequence {@code i}, as BLAST prints it as {@code sseqid}. */
	public String id(int i) {
		return ids.get(i);
	}

	/** Index in {@link #residues()} of the first residue of sequence {@code i}. */
	public int start(int i) {
		return starts[i];
	}

	/** Length of sequence {@code i}. */
	public int length(int i) {
		return starts[i + 1] - starts[i] - 1;
	}

	/** All residues in the {@code .psq} layout; not a copy. */
	public byte[] residues() {
		return residues;
	}

	/** A copy of sequence {@code i} as NCBIstdaa codes. */
	public byte[] sequence(int i) {
		return java.util.Arrays.copyOfRange(residues, starts[i], starts[i] + length(i));
	}

	/** Sequence {@code i} as letters. */
	public String letters(int i) {
		StringBuilder sb = new StringBuilder(length(i));
		for (int k = starts[i], end = starts[i] + length(i); k < end; k++) {
			sb.append(Alphabet.decode(residues[k]));
		}
		return sb.toString();
	}
}
