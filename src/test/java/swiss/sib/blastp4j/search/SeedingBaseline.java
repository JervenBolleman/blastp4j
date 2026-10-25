package swiss.sib.blastp4j.search;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
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
import swiss.sib.blastp4j.seed.SeedParameters;
import swiss.sib.blastp4j.seed.UngappedHsp;
import swiss.sib.blastp4j.stats.ScoringParameters;
import swiss.sib.blastp4j.stats.Statistics;

/**
 * Step 5: the seeding stage for the 24 reference queries against Swiss-Prot, compared with the blastp
 * {@code -comp_based_stats 0} reference.
 * <p>
 * BLAST starts each gapped extension at a point inside an ungapped alignment of the seeding stage, and the final
 * alignment goes through that point. So every reported BLAST alignment must cross at least one of the ungapped
 * alignments found here (same query and subject position on the same diagonal), and every reported subject must reach
 * the gapped stage. The run checks both, and, when {@code target/exhaustive.tsv} from {@link ExhaustiveBaseline} is
 * present, how many of the Smith-Waterman hits that BLAST does not report have no seed (for the queries below BLAST's
 * 500-subject limit).
 * <p>
 * <p>
 * blastp {@code -ungapped} reports the ungapped extensions themselves. Its cutoff is at most the one of the gapped
 * search, so each of its alignments that scores at least our cutoff must be exactly one of the seeds found here: same
 * start, length and score ({@code exact} column).
 * <p>
 * Run from the project directory with the test classpath; the optional argument is the number of threads.
 */
public final class SeedingBaseline {
	public static void main(String[] args) throws Exception {
		int threads = args.length > 0 ? Integer.parseInt(args[0]) : Runtime.getRuntime().availableProcessors();
		Path dir = Path.of(System.getProperty("blastp4j.sprot", System.getProperty("user.home") + "/data/uniprot"));
		SequenceDatabase db = FastaReader.read(dir.resolve("uniprot_sprot.fasta"));
		SequenceDatabase queries = FastaReader.read(Path.of("src/test/resources/queries/uniprot_human.fasta"));
		int minLength = Integer.MAX_VALUE;
		Map<String, Integer> subjectIndex = new HashMap<>();
		for (int i = 0; i < db.size(); i++) {
			minLength = Math.min(minLength, db.length(i));
			subjectIndex.put(db.id(i), i);
		}

		// All BLAST HSPs per query and subject.
		Map<String, Map<String, List<String[]>>> blast = new HashMap<>();
		for (String[] r : ReferenceData.load("human_vs_sprot.cbs0.tsv.gz")) {
			blast.computeIfAbsent(r[ReferenceData.QSEQID], q -> new HashMap<>())
					.computeIfAbsent(r[ReferenceData.SSEQID], s -> new ArrayList<>()).add(r);
		}
		// Ungapped BLAST alignments per query: "subject qstart sstart length score", 0-based.
		Map<String, List<String[]>> ungapped = new HashMap<>();
		for (String[] r : ReferenceData.load("human_vs_sprot.ungapped.tsv.gz")) {
			ungapped.computeIfAbsent(r[ReferenceData.QSEQID], q -> new ArrayList<>()).add(r);
		}
		// Smith-Waterman hits per query, when the exhaustive baseline has been run.
		Map<String, Set<String>> exhaustive = new HashMap<>();
		Path sw = Path.of("target/exhaustive.tsv");
		if (Files.exists(sw)) {
			try (BufferedReader in = Files.newBufferedReader(sw)) {
				for (String line; (line = in.readLine()) != null;) {
					String[] f = line.split("\t");
					exhaustive.computeIfAbsent(f[0], q -> new HashSet<>()).add(f[1]);
				}
			}
		}

		Statistics statistics = new Statistics(ScoringParameters.BLOSUM62_11_1, db.totalResidues());
		SeedSearch search = new SeedSearch(db);
		System.out.printf("%-22s %5s %5s %6s %6s %6s %7s %6s %7s %6s %6s %6s %6s %5s%n", "query", "xdrop", "trig",
				"cutoff", "seeded", "blast", "covered", "hsps", "crossed", "swOnly", "noSeed", "ungap", "exact", "sec");
		long totalNanos = 0;
		int[] totals = new int[9];
		List<String> problems = new ArrayList<>();
		try (ExecutorService pool = Executors.newFixedThreadPool(threads)) {
			for (int q = 0; q < queries.size(); q++) {
				String qid = queries.id(q);
				byte[] query = queries.sequence(q);
				SeedParameters p = SeedParameters.blastp(query, statistics, 10, minLength);
				long t0 = System.nanoTime();
				List<SeedSearch.Seeds> seeds = search.search(query, p, pool, threads);
				long nanos = System.nanoTime() - t0;
				totalNanos += nanos;

				Map<Integer, List<UngappedHsp>> bySubject = new HashMap<>();
				for (SeedSearch.Seeds s : seeds) {
					bySubject.put(s.subject(), s.hsps());
				}
				Map<String, List<String[]>> reference = blast.getOrDefault(qid, Map.of());
				int covered = 0, hsps = 0, crossed = 0;
				for (Map.Entry<String, List<String[]>> e : reference.entrySet()) {
					List<UngappedHsp> ours = bySubject.get(subjectIndex.get(e.getKey()));
					if (ours != null) {
						covered++;
					} else {
						problems.add(qid + " " + e.getKey() + ": no seed");
					}
					for (String[] r : e.getValue()) {
						hsps++;
						if (ours != null && crosses(ours, r)) {
							crossed++;
						} else if (ours != null) {
							problems.add(qid + " " + e.getKey() + " " + r[6] + "-" + r[7] + "/" + r[8] + "-" + r[9]
									+ ": no seed crosses the alignment");
						}
					}
				}
				// Smith-Waterman hits BLAST does not report, below the 500-subject limit only.
				int swOnly = 0, noSeed = 0;
				if (reference.size() < 500) {
					for (String sid : exhaustive.getOrDefault(qid, Set.of())) {
						if (!reference.containsKey(sid)) {
							swOnly++;
							if (!bySubject.containsKey(subjectIndex.get(sid))) {
								noSeed++;
							}
						}
					}
				}
				// Ungapped BLAST alignments at our cutoff, and how many are exactly one of our seeds.
				Set<String> ourKeys = new HashSet<>();
				for (SeedSearch.Seeds s : seeds) {
					for (UngappedHsp h : s.hsps()) {
						ourKeys.add(db.id(s.subject()) + " " + h.queryStart() + " " + h.subjectStart() + " "
								+ h.length() + " " + h.score());
					}
				}
				int ungap = 0, exact = 0;
				for (String[] r : ungapped.getOrDefault(qid, List.of())) {
					int score = Integer.parseInt(r[ReferenceData.SCORE]);
					if (score < p.cutoff()) {
						continue;
					}
					ungap++;
					int qs = Integer.parseInt(r[6]) - 1, ss = Integer.parseInt(r[8]) - 1;
					int length = Integer.parseInt(r[7]) - qs;
					if (ourKeys.contains(r[1] + " " + qs + " " + ss + " " + length + " " + score)) {
						exact++;
					} else {
						problems.add(qid + " " + r[1] + " ungapped " + r[6] + "-" + r[7] + "/" + r[8] + "-" + r[9]
								+ " score " + score + ": not among our seeds");
					}
				}
				System.out.printf("%-22s %5d %5d %6d %6d %6d %7d %6d %7d %6d %6d %6d %6d %5.1f%n", qid, p.xDrop(),
						p.gapTrigger(), p.cutoff(), seeds.size(), reference.size(), covered, hsps, crossed, swOnly,
						noSeed, ungap, exact, nanos / 1e9);
				int[] row = { seeds.size(), reference.size(), covered, hsps, crossed, swOnly, noSeed, ungap, exact };
				for (int i = 0; i < row.length; i++) {
					totals[i] += row[i];
				}
			}
		}
		System.out.printf("%-22s %5s %5s %6s %6d %6d %7d %6d %7d %6d %6d %6d %6d %5.1f%n", "total", "", "", "",
				totals[0], totals[1], totals[2], totals[3], totals[4], totals[5], totals[6], totals[7], totals[8],
				totalNanos / 1e9);
		System.out.println(problems.size() + " problems");
		problems.stream().limit(30).forEach(s -> System.out.println("  " + s));
	}

	/** Whether an ungapped alignment shares a (query, subject) position with the BLAST alignment's diagonal span. */
	static boolean crosses(List<UngappedHsp> hsps, String[] r) {
		int qs = Integer.parseInt(r[6]) - 1, qe = Integer.parseInt(r[7]) - 1;
		int ss = Integer.parseInt(r[8]) - 1, se = Integer.parseInt(r[9]) - 1;
		for (UngappedHsp h : hsps) {
			int from = Math.max(0, Math.max(qs - h.queryStart(), ss - h.subjectStart()));
			int to = Math.min(h.length() - 1, Math.min(qe - h.queryStart(), se - h.subjectStart()));
			if (from <= to) {
				return true;
			}
		}
		return false;
	}
}
