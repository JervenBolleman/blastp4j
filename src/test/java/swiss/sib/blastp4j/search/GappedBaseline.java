package swiss.sib.blastp4j.search;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import swiss.sib.blastp4j.ReferenceData;
import swiss.sib.blastp4j.align.Hsp;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.output.TabularFormat;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Steps 6 and 7: the full blastp search for the 24 reference queries against Swiss-Prot, compared line by line with
 * the blastp reference (all 16 columns, including the alignment as BTOP): without composition-based statistics
 * ({@code cbs0}, the default) or with blastp's default mode 2 ({@code cbs2}). Writes our output to
 * {@code target/blastp4j.<mode>.tsv}. Arguments: [threads] [cbs0|cbs2].
 * <p>
 * Per query: lines in both, lines only in BLAST's output, lines only in ours, and whether the lines in both are in the
 * same order. Run from the project directory with the test classpath; the optional argument is the number of threads.
 */
public final class GappedBaseline {
	public static void main(String[] args) throws Exception {
		int threads = args.length > 0 ? Integer.parseInt(args[0]) : Runtime.getRuntime().availableProcessors();
		String mode = args.length > 1 ? args[1] : "cbs0";
		boolean compositionBased = mode.equals("cbs2");
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		SequenceDatabase db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		Map<String, List<String>> reference = new HashMap<>();
		for (String[] r : ReferenceData.load("human_vs_sprot." + mode + ".tsv.gz")) {
			reference.computeIfAbsent(r[0], q -> new ArrayList<>()).add(String.join("\t", r));
		}
		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		BlastpSearch search = new BlastpSearch(db, statistics, 10, compositionBased);

		System.out.printf("%-22s %6s %6s %6s %6s %6s %6s%n", "query", "blast", "ours", "same", "blOnly", "ourOnly",
				"sec");
		int[] totals = new int[5];
		long totalNanos = 0;
		List<String> examples = new ArrayList<>();
		try (ExecutorService pool = Executors.newFixedThreadPool(threads);
				PrintWriter out = new PrintWriter(Files.newBufferedWriter(Path.of("target/blastp4j." + mode + ".tsv")))) {
			for (int q = 0; q < queries.size(); q++) {
				String qid = queries.id(q);
				byte[] query = queries.sequence(q);
				long t0 = System.nanoTime();
				List<BlastpSearch.Hit> hits = search.search(query, pool, threads);
				long nanos = System.nanoTime() - t0;
				totalNanos += nanos;
				List<String> ours = new ArrayList<>();
				for (BlastpSearch.Hit h : hits) {
					for (Hsp hsp : h.hsps()) {
						ours.add(TabularFormat.row(qid, query, db, h.subject(), hsp, statistics));
					}
				}
				ours.forEach(out::println);
				List<String> theirs = reference.getOrDefault(qid, List.of());
				List<String> blastOnly = new ArrayList<>(theirs);
				List<String> oursOnly = new ArrayList<>();
				for (String line : ours) {
					if (!blastOnly.remove(line)) {
						oursOnly.add(line);
					}
				}
				int same = ours.size() - oursOnly.size();
				List<String> commonTheirs = new ArrayList<>(theirs);
				commonTheirs.removeAll(blastOnly);
				List<String> commonOurs = new ArrayList<>(ours);
				commonOurs.removeAll(oursOnly);
				boolean order = commonTheirs.equals(commonOurs);
				System.out.printf("%-22s %6d %6d %6d %6d %6d %6.1f%s%n", qid, theirs.size(), ours.size(), same,
						blastOnly.size(), oursOnly.size(), nanos / 1e9, order ? "" : "  (order differs)");
				int[] row = { theirs.size(), ours.size(), same, blastOnly.size(), oursOnly.size() };
				for (int i = 0; i < row.length; i++) {
					totals[i] += row[i];
				}
				for (int i = 0; i < Math.min(2, blastOnly.size()) && examples.size() < 40; i++) {
					examples.add("  blast: " + blastOnly.get(i));
				}
				for (int i = 0; i < Math.min(2, oursOnly.size()) && examples.size() < 40; i++) {
					examples.add("  ours:  " + oursOnly.get(i));
				}
			}
		}
		System.out.printf("%-22s %6d %6d %6d %6d %6d %6.1f%n", "total", totals[0], totals[1], totals[2], totals[3],
				totals[4], totalNanos / 1e9);
		examples.forEach(System.out::println);
	}
}
