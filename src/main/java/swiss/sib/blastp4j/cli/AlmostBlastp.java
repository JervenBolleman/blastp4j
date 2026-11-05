package swiss.sib.blastp4j.cli;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.BlastDbReader;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.output.TabularFormat;
import swiss.sib.blastp4j.search.BlastpSearch;
import swiss.sib.blastp4j.search.QueryBatches;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;
import swiss.sib.blastp4j.util.Arithmetic;

/**
 * The {@code almostblastp} command: {@link Blastp} where the results may differ slightly from blastp's, for changes
 * that make the search faster but not bit for bit BLAST's, such as fused multiply-adds in the floating-point code.
 * {@link Blastp} stays exact. A copy of it that differs only in the search's {@link Arithmetic}: {@link Arithmetic#FAST}
 * (for now in the composition-based stage's target-frequency optimizer).
 *
 * <pre>
 * almostblastp -query queries.fasta -db uniprot_sprot.fasta -outfmt "6 std btop" -num_threads 8 -out hits.tsv
 * </pre>
 *
 * The database is a BLAST protein database made by {@code makeblastdb}, or a FASTA file (plain or gzip-compressed).
 * For a FASTA file the results are those blastp gives for a database built from it with
 * {@code makeblastdb -parse_seqids}. Options of blastp that are not implemented are accepted
 * only with their default value.
 */
public final class AlmostBlastp {
	static final String VERSION = "almostblastp (blastp4j 0.1.0, close to BLAST+ 2.17.0)";

	static final String USAGE = """
			Usage: almostblastp -query FILE -db DB [options]

			Searches protein sequences against a protein database and writes the results as BLAST+ blastp's
			tabular output (-outfmt 6). Results may differ slightly from blastp's for the same parameters (for
			blastp's exact results use blastp4j).

			  -query FILE             query sequences, FASTA (plain or .gz)
			  -db DB                  database: a BLAST protein database (the makeblastdb -out name), or
			                          a FASTA file (plain or .gz)
			  -out FILE               output file (default: standard output)
			  -outfmt "6 [FIELDS]"    tabular output; fields as in blastp: std qseqid sseqid qacc qaccver sacc
			                          saccver qlen slen qstart
			                          qend sstart send qseq sseq evalue bitscore score length pident nident
			                          mismatch positive gapopen gaps ppos qcovhsp btop (default: 6, i.e. std)
			  -evalue X               E-value cutoff (default 10)
			  -max_target_seqs N      number of subjects to report per query (default 500)
			  -max_hsps N             number of HSPs to report per subject (default: all)
			  -comp_based_stats MODE  composition-based statistics: 2, T or D (default, conditional score matrix
			                          adjustment) or 0 or F (none)
			  -num_threads N          threads (default: number of processors)
			  -version                print the version
			  -h, -help               print this help

			Accepted with their default value only: -task blastp, -matrix BLOSUM62, -gapopen 11, -gapextend 1,
			-word_size 3, -threshold 11, -window_size 40, -seg no.
			""";

	/** Options that blastp4j implements only with blastp's default value. */
	static final Map<String, String> FIXED = Map.of("task", "blastp", "matrix", "BLOSUM62", "gapopen", "11",
			"gapextend", "1", "word_size", "3", "threshold", "11", "window_size", "40", "seg", "no");

	/** A usage error. */
	static final class UsageException extends Exception {
		UsageException(String message) {
			super(message);
		}
	}

	/** The parsed command line. */
	record Options(Path query, Path db, Path out, List<TabularFormat.Column> columns, double evalue,
			int maxTargetSeqs, int maxHsps, boolean compositionBased, int threads) {
	}

	static Options parse(String[] args) throws UsageException {
		Map<String, String> values = new HashMap<>();
		for (int i = 0; i < args.length; i++) {
			String a = args[i];
			if (!a.startsWith("-") || a.length() < 2) {
				throw new UsageException("unexpected argument '" + a + "'");
			}
			String name = a.substring(a.startsWith("--") ? 2 : 1);
			if (i + 1 >= args.length) {
				throw new UsageException("option -" + name + " needs a value");
			}
			if (values.put(name, args[++i]) != null) {
				throw new UsageException("option -" + name + " is given twice");
			}
		}
		for (String name : values.keySet()) {
			switch (name) {
			case "query", "db", "out", "outfmt", "evalue", "max_target_seqs", "max_hsps", "comp_based_stats",
					"num_threads" -> {
			}
			default -> {
				String fixed = FIXED.get(name);
				if (fixed == null) {
					throw new UsageException("unknown or unsupported option -" + name);
				}
				if (!fixed.equalsIgnoreCase(values.get(name))) {
					throw new UsageException("-" + name + " " + values.get(name) + " is not supported; only "
							+ fixed + " (the default) is");
				}
			}
			}
		}
		if (!values.containsKey("query") || !values.containsKey("db")) {
			throw new UsageException("-query and -db are required");
		}
		String outfmt = values.getOrDefault("outfmt", "6").trim();
		String[] fmt = outfmt.split("\\s+", 2);
		if (!fmt[0].equals("6")) {
			throw new UsageException("only -outfmt 6 (tabular) is supported, not '" + fmt[0] + "'");
		}
		List<TabularFormat.Column> columns;
		try {
			columns = TabularFormat.parseColumns(fmt.length > 1 ? fmt[1] : "std");
		} catch (IllegalArgumentException e) {
			throw new UsageException(e.getMessage());
		}
		boolean cbs = switch (values.getOrDefault("comp_based_stats", "2")) {
		case "2", "T", "t", "D", "d" -> true;
		case "0", "F", "f" -> false;
		default -> throw new UsageException("-comp_based_stats " + values.get("comp_based_stats")
				+ " is not supported; use 2 (the default) or 0");
		};
		double evalue = number(values, "evalue", 10.0);
		if (!(evalue > 0)) {
			throw new UsageException("-evalue must be positive");
		}
		int maxTargetSeqs = (int) integer(values, "max_target_seqs", 500, 1);
		int maxHsps = (int) integer(values, "max_hsps", 0, 1);
		int threads = (int) integer(values, "num_threads", Runtime.getRuntime().availableProcessors(), 1);
		return new Options(Path.of(values.get("query")), Path.of(values.get("db")),
				values.containsKey("out") ? Path.of(values.get("out")) : null, columns, evalue, maxTargetSeqs, maxHsps,
				cbs, threads);
	}

	private static double number(Map<String, String> values, String name, double dflt) throws UsageException {
		String v = values.get(name);
		if (v == null) {
			return dflt;
		}
		try {
			return Double.parseDouble(v);
		} catch (NumberFormatException e) {
			throw new UsageException("-" + name + " needs a number, not '" + v + "'");
		}
	}

	private static long integer(Map<String, String> values, String name, long dflt, long min) throws UsageException {
		String v = values.get(name);
		if (v == null) {
			return dflt;
		}
		try {
			long n = Long.parseLong(v);
			if (n < min || n > Integer.MAX_VALUE) {
				throw new UsageException("-" + name + " must be at least " + min);
			}
			return n;
		} catch (NumberFormatException e) {
			throw new UsageException("-" + name + " needs a whole number, not '" + v + "'");
		}
	}

	/** Runs the search; returns the exit status. */
	static int run(String[] args, PrintStream stdout, PrintStream stderr) {
		if (args.length == 0 || List.of(args).contains("-h") || List.of(args).contains("-help")
				|| List.of(args).contains("--help")) {
			stdout.print(USAGE);
			return args.length == 0 ? 1 : 0;
		}
		if (List.of(args).contains("-version") || List.of(args).contains("--version")) {
			stdout.println(VERSION);
			return 0;
		}
		Options o;
		try {
			o = parse(args);
		} catch (UsageException e) {
			stderr.println("almostblastp: " + e.getMessage());
			stderr.println("Try 'almostblastp -help'.");
			return 1;
		}
		try {
			search(o, stdout, stderr);
			return 0;
		} catch (IOException e) {
			stderr.println("almostblastp: " + e.getMessage());
			return 2;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			return 3;
		}
	}

	static void search(Options o, PrintStream stdout, PrintStream stderr) throws IOException, InterruptedException {
		if (!Files.isReadable(o.query())) {
			throw new IOException("cannot read query file " + o.query());
		}
		SequenceDatabase queries = FastaReader.read(o.query());
		SequenceDatabase db = openDatabase(o.db());
		try (db; ExecutorService pool = Executors.newFixedThreadPool(o.threads());
				Writer out = o.out() != null ? Files.newBufferedWriter(o.out(), StandardCharsets.US_ASCII)
						: new BufferedWriter(new OutputStreamWriter(stdout, StandardCharsets.US_ASCII))) {
			search(o, queries, db, new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues()), pool, out,
					stderr);
		}
	}

	/** The database at {@code path}: a BLAST protein database, or a FASTA file. */
	static SequenceDatabase openDatabase(Path path) throws IOException {
		if (BlastDbReader.exists(path)) {
			return BlastDbReader.read(path);
		} else if (Files.isReadable(path)) {
			return FastaReader.read(path);
		}
		throw new IOException("cannot read database " + path + " (neither a FASTA file nor a BLAST database)");
	}

	/**
	 * Searches {@code queries} against {@code db} with the options' search parameters, in blastp's query batches, on
	 * {@code o.threads()} tasks of {@code pool}; writes the tabular rows to {@code out} and the warnings about skipped
	 * queries to {@code stderr}.
	 */
	static void search(Options o, SequenceDatabase queries, SequenceDatabase db, Statistics statistics,
			ExecutorService pool, Writer out, PrintStream stderr) throws IOException, InterruptedException {
		// Where almostblastp differs from blastp4j: the floating-point code may compute in its faster way.
		BlastpSearch search = new BlastpSearch(db, statistics, o.evalue(), o.compositionBased(), o.maxTargetSeqs(),
				o.maxHsps(), Arithmetic.FAST);
		int[] lengths = new int[queries.size()];
		for (int q = 0; q < lengths.length; q++) {
			lengths[q] = queries.length(q);
		}
		QueryBatches batches = QueryBatches.blastp(lengths);
		for (int first = 0, end; first < queries.size(); first = end) {
			// BLAST's batch: queries until it holds the batch size in residues. The preliminary stage scans the
			// database once for the batch; the later stages go query by query.
			long residues = 0;
			for (end = first; end < queries.size() && residues < batches.batchSize(); end++) {
				residues += lengths[end];
			}
			// As blastp, skip a query whose composition gives no Karlin-Altschul parameters.
			String[] skipped = new String[end - first];
			List<byte[]> batch = new ArrayList<>();
			for (int q = first; q < end; q++) {
				try {
					search.seedParameters(queries.sequence(q));
					search.cutoff(queries.sequence(q));
					batch.add(queries.sequence(q));
				} catch (IllegalArgumentException e) {
					skipped[q - first] = e.getMessage();
				}
			}
			List<List<BlastpSearch.Hit>> prelim = batch.size() > 1 ? search.prelim(batch, pool, o.threads()) : null;
			int inBatch = 0;
			for (int q = first; q < end; q++) {
				byte[] query = queries.sequence(q);
				boolean earlyTermination = batches.next();
				if (skipped[q - first] != null) {
					batches.done(false);
					stderr.println("almostblastp: warning: query " + queries.id(q) + " skipped: " + skipped[q - first]);
					continue;
				}
				List<BlastpSearch.Hit> hits;
				try {
					BlastpSearch.Result result = prelim != null
							? search.finish(query, prelim.get(inBatch++), pool, o.threads(), earlyTermination)
							: search.search(query, pool, o.threads(), earlyTermination);
					batches.done(result.filledToCutoff());
					hits = result.hits();
				} catch (IllegalArgumentException e) {
					batches.done(false);
					stderr.println("almostblastp: warning: query " + queries.id(q) + " skipped: " + e.getMessage());
					continue;
				}
				for (BlastpSearch.Hit h : hits) {
					for (Hsp hsp : h.hsps()) {
						out.write(TabularFormat.row(o.columns(), queries.name(q), query, db, h.subject(), hsp,
								statistics));
						out.write('\n');
					}
				}
				out.flush();
			}
		}
	}

	public static void main(String[] args) {
		Locale.setDefault(Locale.ROOT);
		System.exit(run(args, System.out, System.err));
	}
}
