package swiss.sib.blastp4j.db;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Reads a BLAST protein database, as made by {@code makeblastdb -dbtype prot} (format version 4 or 5), into a
 * {@link SequenceDatabase}.
 * <p>
 * A volume consists of:
 * <ul>
 * <li>{@code .pin}, the index: version, type, title, date, number of sequences, total and maximum length, then the
 * offsets of every header in {@code .phr} and of every sequence in {@code .psq}. Integers are big-endian, except the
 * total length (64 bits, little-endian);</li>
 * <li>{@code .psq}, the sequences as NCBIstdaa codes separated by a 0 byte: the in-memory layout of
 * {@link SequenceDatabase}, which is loaded as is;</li>
 * <li>{@code .phr}, the headers: one {@code Blast-def-line-set} per sequence in binary ASN.1, decoded only when an
 * identifier or title is asked for.</li>
 * </ul>
 * A multi-volume database has a {@code .pal} alias file whose {@code DBLIST} names the volumes; alias files that
 * restrict the database ({@code OIDLIST}, {@code GILIST}, {@code SEQIDLIST}, {@code TAXIDLIST}...) are not
 * supported. The other files (LMDB, taxonomy, masks) are not needed.
 */
public final class BlastDbReader {
	private BlastDbReader() {
	}

	/** Whether {@code prefix} names a BLAST protein database ({@code prefix.pin} or {@code prefix.pal} exists). */
	public static boolean exists(Path prefix) {
		return Files.exists(sibling(prefix, ".pin")) || Files.exists(sibling(prefix, ".pal"));
	}

	private static Path sibling(Path prefix, String extension) {
		return prefix.resolveSibling(prefix.getFileName() + extension);
	}

	/** Reads the database {@code prefix} (the name given to {@code makeblastdb -out}, without extension). */
	public static SequenceDatabase read(Path prefix) throws IOException {
		List<Path> volumes = new ArrayList<>();
		collectVolumes(prefix, volumes, 0);
		List<Volume> read = new ArrayList<>();
		long total = 1;
		for (Path v : volumes) {
			Volume vol = readVolume(v);
			read.add(vol);
			total += vol.sequences.length - 1;
		}
		if (total > Integer.MAX_VALUE - 16) {
			throw new IOException("database too large for one array: " + total + " bytes of sequence data");
		}
		if (read.size() == 1) {
			Volume v = read.get(0);
			return new SequenceDatabase(v.sequences, v.seqOffsets, new DbHeaders(List.of(v)));
		}
		// Concatenate the volumes; each .psq starts with the sentinel that ends the previous one.
		byte[] residues = new byte[(int) total];
		int count = 0;
		for (Volume v : read) {
			count += v.count;
		}
		int[] starts = new int[count + 1];
		int pos = 1, k = 0;
		for (Volume v : read) {
			int base = pos - 1;
			System.arraycopy(v.sequences, 1, residues, pos, v.sequences.length - 1);
			for (int i = 0; i < v.count; i++) {
				starts[k++] = v.seqOffsets[i] + base;
			}
			pos += v.sequences.length - 1;
		}
		starts[count] = pos;
		return new SequenceDatabase(residues, starts, new DbHeaders(read));
	}

	/** Follows alias files ({@code .pal}) down to volumes ({@code .pin}). */
	private static void collectVolumes(Path prefix, List<Path> volumes, int depth) throws IOException {
		if (depth > 16) {
			throw new IOException("alias files nested too deeply at " + prefix);
		}
		Path pal = sibling(prefix, ".pal");
		if (Files.exists(pal)) {
			List<String> dblist = null;
			for (String line : Files.readAllLines(pal, StandardCharsets.ISO_8859_1)) {
				String t = line.trim();
				if (t.isEmpty() || t.startsWith("#")) {
					continue;
				}
				String[] kv = t.split("\\s+", 2);
				String key = kv[0].toUpperCase(Locale.ROOT);
				switch (key) {
				case "DBLIST" -> dblist = splitList(kv.length > 1 ? kv[1] : "");
				case "TITLE", "NSEQ", "LENGTH", "MAXOID", "MEMB_BIT", "DATE" -> {
				}
				default -> throw new IOException(pal + ": alias key " + key + " is not supported");
				}
			}
			if (dblist == null || dblist.isEmpty()) {
				throw new IOException(pal + ": no DBLIST");
			}
			for (String name : dblist) {
				Path p = Path.of(name);
				collectVolumes(p.isAbsolute() ? p : prefix.resolveSibling(name), volumes, depth + 1);
			}
			return;
		}
		if (!Files.exists(sibling(prefix, ".pin"))) {
			throw new IOException("no BLAST protein database " + prefix + " (" + prefix.getFileName()
					+ ".pin or .pal not found)");
		}
		volumes.add(prefix);
	}

	/** Splits a DBLIST value: names separated by spaces, possibly in double quotes. */
	private static List<String> splitList(String s) {
		List<String> out = new ArrayList<>();
		StringBuilder cur = new StringBuilder();
		boolean quoted = false;
		for (char c : s.toCharArray()) {
			if (c == '"') {
				quoted = !quoted;
			} else if (Character.isWhitespace(c) && !quoted) {
				if (!cur.isEmpty()) {
					out.add(cur.toString());
					cur.setLength(0);
				}
			} else {
				cur.append(c);
			}
		}
		if (!cur.isEmpty()) {
			out.add(cur.toString());
		}
		return out;
	}

	/** One volume: its sequence file, offsets and headers. */
	static final class Volume {
		final String title;
		final int count;
		final long totalLength;
		final byte[] sequences;
		final int[] seqOffsets;
		final int[] headerOffsets;
		final byte[] headers;

		Volume(String title, int count, long totalLength, byte[] sequences, int[] seqOffsets, int[] headerOffsets,
				byte[] headers) {
			this.title = title;
			this.count = count;
			this.totalLength = totalLength;
			this.sequences = sequences;
			this.seqOffsets = seqOffsets;
			this.headerOffsets = headerOffsets;
			this.headers = headers;
		}
	}

	static Volume readVolume(Path prefix) throws IOException {
		Path pin = sibling(prefix, ".pin");
		ByteBuffer in = ByteBuffer.wrap(Files.readAllBytes(pin)).order(ByteOrder.BIG_ENDIAN);
		int version = in.getInt();
		if (version != 4 && version != 5) {
			throw new IOException(pin + ": unsupported format version " + version);
		}
		int type = in.getInt();
		if (type != 1) {
			throw new IOException(pin + ": not a protein database");
		}
		if (version == 5) {
			in.getInt(); // volume number
		}
		String title = string(in);
		if (version == 5) {
			string(in); // LMDB file name
		}
		string(in); // date
		int count = in.getInt();
		long totalLength = in.order(ByteOrder.LITTLE_ENDIAN).getLong();
		in.order(ByteOrder.BIG_ENDIAN);
		in.getInt(); // maximum length
		int[] headerOffsets = new int[count + 1];
		for (int i = 0; i <= count; i++) {
			headerOffsets[i] = in.getInt();
		}
		int[] seqOffsets = new int[count + 1];
		for (int i = 0; i <= count; i++) {
			seqOffsets[i] = in.getInt();
		}
		byte[] sequences = Files.readAllBytes(sibling(prefix, ".psq"));
		if (seqOffsets[count] != sequences.length) {
			throw new IOException(prefix + ".psq: " + sequences.length + " bytes, the index says " + seqOffsets[count]);
		}
		byte[] headers = Files.readAllBytes(sibling(prefix, ".phr"));
		if (headerOffsets[count] != headers.length) {
			throw new IOException(prefix + ".phr: " + headers.length + " bytes, the index says " + headerOffsets[count]);
		}
		return new Volume(title, count, totalLength, sequences, seqOffsets, headerOffsets, headers);
	}

	/** A length-prefixed string; makeblastdb pads some with zero bytes for alignment. */
	private static String string(ByteBuffer in) {
		int n = in.getInt();
		byte[] b = new byte[n];
		in.get(b);
		int end = n;
		while (end > 0 && b[end - 1] == 0) {
			end--;
		}
		return new String(b, 0, end, StandardCharsets.ISO_8859_1);
	}

	/** Headers decoded from the {@code .phr} data of the volumes, on demand. */
	static final class DbHeaders implements Headers {
		private final List<Volume> volumes;
		private final int[] firstOid;
		private final int size;

		DbHeaders(List<Volume> volumes) {
			this.volumes = volumes;
			firstOid = new int[volumes.size()];
			int n = 0;
			for (int v = 0; v < volumes.size(); v++) {
				firstOid[v] = n;
				n += volumes.get(v).count;
			}
			size = n;
		}

		@Override
		public int size() {
			return size;
		}

		/** The identifiers and title of the first definition line of sequence {@code i}. */
		private record Defline(List<SeqId> ids, String title) {
		}

		private Defline defline(int i) {
			int v = volumes.size() - 1;
			while (firstOid[v] > i) {
				v--;
			}
			Volume vol = volumes.get(v);
			int oid = i - firstOid[v];
			ByteBuffer b = ByteBuffer.wrap(vol.headers, vol.headerOffsets[oid],
					vol.headerOffsets[oid + 1] - vol.headerOffsets[oid]);
			// Blast-def-line-set ::= SEQUENCE OF Blast-def-line; the first one is the one BLAST prints.
			Asn1.Node first = Asn1.read(b).children().get(0);
			// Blast-def-line ::= SEQUENCE { title [0], seqid [1] SEQUENCE OF Seq-id, taxid [2], ... }
			Asn1.Node t = first.child(0);
			String title = t == null ? "" : t.inner().string();
			List<SeqId> ids = new ArrayList<>();
			Asn1.Node seqids = first.child(1);
			if (seqids != null) {
				for (Asn1.Node choice : seqids.inner().children()) {
					ids.add(SeqId.fromAsn(choice));
				}
			}
			return new Defline(ids, title);
		}

		/**
		 * Whether the database was made without {@code -parse_seqids}: its only identifier is
		 * {@code gnl|BL_ORD_ID|n}, and BLAST prints the first word of the title instead.
		 */
		private static boolean ordinalOnly(List<SeqId> ids) {
			return ids.size() == 1 && ids.get(0).kind() == SeqId.Kind.GENERAL
					&& ids.get(0).fasta().startsWith("gnl|BL_ORD_ID|");
		}

		private static String firstWord(String s) {
			int end = 0;
			while (end < s.length() && s.charAt(end) > ' ') {
				end++;
			}
			return s.substring(0, end);
		}

		@Override
		public String name(int i) {
			return seqid(i);
		}

		@Override
		public String seqid(int i) {
			Defline d = defline(i);
			return ordinalOnly(d.ids) || d.ids.isEmpty() ? firstWord(d.title) : SeqId.seqid(d.ids);
		}

		@Override
		public String accession(int i, boolean withVersion) {
			Defline d = defline(i);
			return ordinalOnly(d.ids) || d.ids.isEmpty() ? firstWord(d.title) : SeqId.accession(d.ids, withVersion);
		}

		@Override
		public String title(int i) {
			return defline(i).title;
		}
	}
}
