package swiss.sib.blastp4j.search;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import swiss.sib.blastp4j.ReferenceData;
import swiss.sib.blastp4j.db.FastaReader;
import swiss.sib.blastp4j.db.SequenceDatabase;
import swiss.sib.blastp4j.output.TabularFormat;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Step 4: the exhaustive Smith-Waterman baseline for the 24 reference queries against Swiss-Prot, compared with the
 * blastp {@code -comp_based_stats 0} reference. Prints one line per query and writes all hits to
 * {@code target/exhaustive.tsv}.
 * <p>
 * Run from the project directory with the test classpath; the optional argument is the number of threads.
 * <p>
 * Per query: BLAST subjects (with the best HSP per subject), how many the exhaustive search finds and how its score
 * compares, and the exhaustive hits BLAST does not report. BLAST keeps at most 500 subjects, so for a query that
 * reaches the cap only the exhaustive hits ranked within the first 500 count as missed by BLAST.
 */
public final class ExhaustiveBaseline {
	static final int MAX_TARGET_SEQS = 500;

	public static void main(String[] args) throws Exception {
		int threads = args.length > 0 ? Integer.parseInt(args[0]) : Runtime.getRuntime().availableProcessors();
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		long t0 = System.nanoTime();
		SequenceDatabase db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		System.out.printf("loaded %d subjects in %.1f s%n", db.size(), (System.nanoTime() - t0) / 1e9);

		// Best BLAST HSP per query and subject.
		Map<String, Map<String, String[]>> blast = new HashMap<>();
		for (String[] r : ReferenceData.load("human_vs_sprot.cbs0.tsv.gz")) {
			blast.computeIfAbsent(r[ReferenceData.QSEQID], q -> new HashMap<>()).merge(r[ReferenceData.SSEQID], r,
					(a, b) -> Integer.parseInt(a[ReferenceData.SCORE]) >= Integer.parseInt(b[ReferenceData.SCORE]) ? a
							: b);
		}

		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		ExhaustiveSearch search = new ExhaustiveSearch(db, statistics, 10);
		long cells = 0, totalNanos = 0;
		int sumBlast = 0, sumFound = 0, sumHigher = 0, sumSameEnds = 0, sumHits = 0, sumExtra = 0;
		System.out.printf("%-22s %5s %6s %6s %6s %6s %6s %6s %5s %8s%n", "query", "qlen", "blast", "found", "higher",
				"ends=", "sw", "extra", "cap", "seconds");
		try (ExecutorService pool = Executors.newFixedThreadPool(threads);
				PrintWriter out = new PrintWriter(Files.newBufferedWriter(Path.of("target/exhaustive.tsv")))) {
			for (int q = 0; q < queries.size(); q++) {
				String qid = queries.id(q);
				long start = System.nanoTime();
				List<ExhaustiveSearch.Hit> hits = search.search(queries.sequence(q), pool, threads);
				long nanos = System.nanoTime() - start;
				totalNanos += nanos;
				cells += (long) queries.length(q) * db.totalResidues();

				Map<String, String[]> reference = blast.getOrDefault(qid, Map.of());
				boolean capped = reference.size() >= MAX_TARGET_SEQS;
				Set<String> seen = new HashSet<>();
				int found = 0, higher = 0, sameEnds = 0, extra = 0, rank = 0;
				for (ExhaustiveSearch.Hit h : hits) {
					String sid = db.id(h.subject());
					seen.add(sid);
					out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d%n", qid, sid, h.queryStart() + 1,
							h.queryEnd() + 1, h.subjectStart() + 1, h.subjectEnd() + 1,
							TabularFormat.evalue(h.evalue()), TabularFormat.bitScore(statistics.bitScore(h.score())),
							h.score(), queries.length(q), db.length(h.subject()));
					String[] r = reference.get(sid);
					rank++;
					if (r == null) {
						if (!capped || rank <= MAX_TARGET_SEQS) {
							extra++;
						}
						continue;
					}
					found++;
					int score = Integer.parseInt(r[ReferenceData.SCORE]);
					if (h.score() < score) {
						throw new AssertionError(qid + " " + sid + ": exhaustive " + h.score() + " < BLAST " + score);
					}
					if (h.score() > score) {
						higher++;
					} else if (h.queryEnd() + 1 == Integer.parseInt(r[7]) && h.subjectEnd() + 1 == Integer.parseInt(r[9])
							&& h.queryStart() + 1 == Integer.parseInt(r[6]) && h.subjectStart() + 1 == Integer.parseInt(r[8])) {
						sameEnds++;
					}
				}
				for (String sid : reference.keySet()) {
					if (!seen.contains(sid)) {
						System.out.println("  not found: " + sid);
					}
				}
				System.out.printf("%-22s %5d %6d %6d %6d %6d %6d %6d %5s %8.1f%n", qid, queries.length(q),
						reference.size(), found, higher, sameEnds, hits.size(), extra, capped ? "yes" : "", nanos / 1e9);
				sumBlast += reference.size();
				sumFound += found;
				sumHigher += higher;
				sumSameEnds += sameEnds;
				sumHits += hits.size();
				sumExtra += extra;
			}
		}
		System.out.printf("%-22s %5s %6d %6d %6d %6d %6d %6d %5s %8.1f%n", "total", "", sumBlast, sumFound, sumHigher,
				sumSameEnds, sumHits, sumExtra, "", totalNanos / 1e9);
		System.out.printf("%d threads, %.2f GCUPS%n", threads, cells / (double) totalNanos);
	}
}
