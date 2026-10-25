package swiss.sib.blastp4j.stats;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.ReferenceData;
import swiss.sib.blastp4j.output.TabularFormat;

/**
 * Recomputes the bit score and E-value columns of the BLAST+ reference (without composition-based statistics) from
 * the raw score, query length and subject length of each hit, and compares the formatted values exactly.
 */
class StatisticsReferenceTest {

	@Test
	void reproducesBlastpBitScoresAndEvalues() throws IOException {
		Statistics stats = new Statistics(ScoringParameters.BLOSUM62_11_1, ReferenceData.SPROT_RESIDUES);
		List<String[]> rows = ReferenceData.load("human_vs_sprot.cbs0.tsv.gz");
		List<String> wrong = new ArrayList<>();
		for (String[] r : rows) {
			int score = Integer.parseInt(r[ReferenceData.SCORE]);
			int qlen = Integer.parseInt(r[ReferenceData.QLEN]);
			int slen = Integer.parseInt(r[ReferenceData.SLEN]);
			String evalue = TabularFormat.evalue(stats.evalue(score, qlen, slen));
			String bits = TabularFormat.bitScore(stats.bitScore(score));
			if (!evalue.equals(r[ReferenceData.EVALUE]) || !bits.equals(r[ReferenceData.BITSCORE])) {
				wrong.add(r[ReferenceData.QSEQID] + " " + r[ReferenceData.SSEQID] + " score " + score + ": expected "
						+ r[ReferenceData.EVALUE] + " / " + r[ReferenceData.BITSCORE] + ", got " + evalue + " / " + bits);
			}
		}
		assertEquals(15_107, rows.size());
		assertEquals(List.of(), wrong.subList(0, Math.min(10, wrong.size())),
				wrong.size() + " of " + rows.size() + " hits differ");
	}
}
