package swiss.sib.blastp4j.compo;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;

import swiss.sib.blastp4j.stats.UngappedKarlin;

/**
 * The adjusted matrices are those of BLAST+'s own composition adjustment library: 306 query/subject compositions
 * (three queries against 100 random Swiss-Prot subjects each, SEG-masked or not, plus UBB_BOVIN), run through
 * {@code Blast_AdjustScores} compiled from the BLAST+ sources by {@code reference/compo/build_harness.sh}.
 */
class CompositionAdjustmentTest {
	static List<String> lines(String name) throws IOException {
		try (BufferedReader in = new BufferedReader(new InputStreamReader(
				new GZIPInputStream(CompositionAdjustmentTest.class.getResourceAsStream("/compo/" + name)),
				StandardCharsets.US_ASCII))) {
			return in.lines().toList();
		}
	}

	static double[] doubles(String[] f, int n) {
		double[] v = new double[n];
		for (int i = 0; i < n; i++) {
			v[i] = Double.parseDouble(f[i]);
		}
		return v;
	}

	@Test
	void sameMatricesAsBlast() throws IOException {
		List<String> cases = lines("adjust_cases.txt.gz"), blast = lines("adjust_blast.txt.gz");
		int n = Integer.parseInt(cases.get(0).trim());
		double lambda = Double.parseDouble(cases.get(1).trim());
		CompositionAdjustment.MatrixInfo info = new CompositionAdjustment.MatrixInfo(UngappedKarlin.ideal().lambda(),
				32);
		assertEquals(lambda, info.ungappedLambda, "ideal lambda / 32");
		CompositionAdjustment.Workspace ws = new CompositionAdjustment.Workspace();
		int[] rules = new int[6];
		for (int c = 0; c < n; c++) {
			String[] q = cases.get(2 + 2 * c).trim().split(" "), s = cases.get(3 + 2 * c).trim().split(" ");
			var qc = new CompositionAdjustment.Composition(doubles(q, 28), Integer.parseInt(q[28]));
			var sc = new CompositionAdjustment.Composition(doubles(s, 28), Integer.parseInt(s[28]));
			int[][] m = new int[28][28];
			var adj = CompositionAdjustment.adjustScores(m, qc, Integer.parseInt(q[29]), sc, Integer.parseInt(s[29]),
					info, ws);

			String[] header = blast.get(29 * c).trim().split(" ");
			assertEquals(Integer.parseInt(header[0]), adj.status(), "case " + c + " status");
			assertEquals(Integer.parseInt(header[1]), adj.rule().ordinal() - 1, "case " + c + " rule");
			double ratio = Double.parseDouble(header[2]);
			assertEquals(ratio, adj.lambdaRatio(), 1e-12 * Math.abs(ratio), "case " + c + " lambda ratio");
			for (int i = 0; i < 28; i++) {
				int[] expected = Arrays.stream(blast.get(29 * c + 1 + i).trim().split(" "))
						.mapToInt(Integer::parseInt).toArray();
				assertEquals(Arrays.toString(expected), Arrays.toString(m[i]), "case " + c + " row " + i);
			}
			rules[adj.rule().ordinal()]++;
		}
		// Both kinds of adjustment occur in the data.
		List<String> used = new ArrayList<>();
		for (CompositionAdjustment.Rule r : CompositionAdjustment.Rule.values()) {
			if (rules[r.ordinal()] > 0) {
				used.add(r + "=" + rules[r.ordinal()]);
			}
		}
		assertEquals(2, used.size(), used.toString());
	}
}
