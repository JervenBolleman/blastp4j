package swiss.sib.blastp4j.db;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;

/**
 * Reads protein FASTA (plain or gzip-compressed) into a {@link SequenceDatabase}, encoding residues as NCBIstdaa the
 * way {@code makeblastdb} does for the letters it accepts. The id of a sequence is the first word of its header.
 * Parses bytes directly: whitespace inside sequences is skipped, letters outside the alphabet become X.
 */
public final class FastaReader {
	private static final int BUFFER = 1 << 16;

	private FastaReader() {
	}

	public static SequenceDatabase read(Path file) throws IOException {
		long size = Files.size(file);
		boolean gzip = file.getFileName().toString().endsWith(".gz");
		try (InputStream raw = Files.newInputStream(file);
				InputStream in = gzip ? new GZIPInputStream(raw, BUFFER) : new BufferedInputStream(raw, BUFFER)) {
			// A plain file has about as many residues as bytes; gzip-compressed protein FASTA expands about 3 times.
			return read(in, gzip ? size * 3 : size);
		}
	}

	public static SequenceDatabase read(InputStream in, long sizeHint) throws IOException {
		byte[] residues = new byte[(int) Math.min(Integer.MAX_VALUE - 16, Math.max(1024, sizeHint + 2))];
		int n = 0;
		residues[n++] = Alphabet.GAP; // leading sentinel
		int[] starts = new int[1024];
		int count = 0;
		FastaHeaders.Builder headers = new FastaHeaders.Builder();
		byte[] header = new byte[256];
		byte[] buf = new byte[BUFFER];
		boolean lineStart = true, inHeader = false, inSequence = false;
		int headerLength = 0;
		for (int read = in.read(buf); read > 0; read = in.read(buf)) {
			for (int k = 0; k < read; k++) {
				byte b = buf[k];
				if (inHeader) {
					if (b == '\n') {
						headers.add(header, 0, headerLength);
						inHeader = false;
						lineStart = true;
					} else {
						if (headerLength == header.length) {
							header = Arrays.copyOf(header, header.length * 2);
						}
						header[headerLength++] = b;
					}
					continue;
				}
				if (b == '\n') {
					lineStart = true;
					continue;
				}
				if (lineStart && b == '>') {
					if (inSequence) {
						residues = ensure(residues, n + 1);
						residues[n++] = Alphabet.GAP; // sentinel after the previous sequence
					}
					if (count + 2 > starts.length) {
						starts = Arrays.copyOf(starts, starts.length * 2);
					}
					starts[count++] = n;
					inHeader = true;
					inSequence = true;
					headerLength = 0;
					lineStart = false;
					continue;
				}
				lineStart = false;
				if (b <= ' ') {
					continue; // \r, spaces, tabs
				}
				if (!inSequence) {
					throw new IOException("sequence data before the first header");
				}
				residues = ensure(residues, n + 1);
				residues[n++] = Alphabet.encode(b);
			}
		}
		if (inHeader) {
			headers.add(header, 0, headerLength); // header without a final newline
		}
		if (inSequence) {
			residues = ensure(residues, n + 1);
			residues[n++] = Alphabet.GAP;
		}
		starts[count] = n;
		return new SequenceDatabase(n == residues.length ? residues : Arrays.copyOf(residues, n),
				Arrays.copyOf(starts, count + 1), headers.build());
	}

	private static byte[] ensure(byte[] array, int size) {
		if (size <= array.length) {
			return array;
		}
		long grown = Math.max(size, (long) array.length + (array.length >> 1));
		if (grown > Integer.MAX_VALUE - 16) {
			throw new IllegalStateException("database too large for one array: " + grown + " bytes");
		}
		return Arrays.copyOf(array, (int) grown);
	}
}
