package swiss.sib.blastp4j.db;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;

/**
 * Headers of a FASTA file, kept as bytes and interpreted as {@code makeblastdb -parse_seqids} does: the first word
 * holds the identifiers ({@link SeqId#parseFasta}) and the rest is the title.
 *
 * @param text    the header lines, one after the other
 * @param offsets header {@code i} is {@code text[offsets[i] .. offsets[i + 1])}
 */
record FastaHeaders(byte[] text, int[] offsets) implements Headers {

	/** Collects header lines (without the '>'). */
	static final class Builder {
		private byte[] text = new byte[1 << 16];
		private int[] offsets = new int[1024];
		private int count, length;

		void add(byte[] line, int from, int to) {
			int n = to - from;
			if (length + n > text.length) {
				text = Arrays.copyOf(text, Math.max(length + n, text.length * 2));
			}
			System.arraycopy(line, from, text, length, n);
			if (count + 2 > offsets.length) {
				offsets = Arrays.copyOf(offsets, offsets.length * 2);
			}
			offsets[count++] = length;
			length += n;
		}

		FastaHeaders build() {
			offsets[count] = length;
			return new FastaHeaders(Arrays.copyOf(text, length), Arrays.copyOf(offsets, count + 1));
		}

		int size() {
			return count;
		}
	}

	@Override
	public int size() {
		return offsets.length - 1;
	}

	private int wordEnd(int i) {
		int end = offsets[i];
		while (end < offsets[i + 1] && text[end] > ' ') {
			end++;
		}
		return end;
	}

	@Override
	public String name(int i) {
		return new String(text, offsets[i], wordEnd(i) - offsets[i], StandardCharsets.ISO_8859_1);
	}

	private List<SeqId> ids(int i) {
		return SeqId.parseFasta(name(i));
	}

	@Override
	public String seqid(int i) {
		return SeqId.seqid(ids(i));
	}

	@Override
	public String accession(int i, boolean withVersion) {
		return SeqId.accession(ids(i), withVersion);
	}

	@Override
	public String title(int i) {
		int from = wordEnd(i), to = offsets[i + 1];
		while (from < to && text[from] <= ' ') {
			from++;
		}
		while (to > from && text[to - 1] <= ' ') {
			to--;
		}
		return new String(text, from, to - from, StandardCharsets.ISO_8859_1);
	}
}
