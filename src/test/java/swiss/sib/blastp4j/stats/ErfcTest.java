package swiss.sib.blastp4j.stats;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

/** erfc against reference values from Python's math.erfc (C library erfc). */
class ErfcTest {
	private static final double[][] VALUES = { //
			{ -6.0, 2.0 }, //
			{ -3.0, 1.9999779095030015 }, //
			{ -1.5, 1.9661051464753108 }, //
			{ -1.0, 1.8427007929497148 }, //
			{ -0.5, 1.5204998778130465 }, //
			{ -0.1, 1.1124629160182848 }, //
			{ 0.0, 1.0 }, //
			{ 1e-08, 0.9999999887162083 }, //
			{ 0.1, 0.887537083981715 }, //
			{ 0.5, 0.4795001221869535 }, //
			{ 1.0, 0.15729920705028513 }, //
			{ 1.2247, 0.08327581479808639 }, //
			{ 1.2248, 0.08325063757190905 }, //
			{ 1.5, 0.033894853524689274 }, //
			{ 2.0, 0.004677734981047266 }, //
			{ 3.0, 2.209049699858544e-05 }, //
			{ 4.0, 1.541725790028002e-08 }, //
			{ 5.0, 1.537459794428035e-12 }, //
			{ 6.0, 2.1519736712498913e-17 }, //
			{ 10.0, 2.088487583762545e-45 }, //
			{ 20.0, 5.395865611607901e-176 }, //
			{ 26.0, 5.663192408856143e-296 } };

	@Test
	void matchesReferenceValues() {
		for (double[] v : VALUES) {
			double expected = v[1];
			assertEquals(expected, Erfc.erfc(v[0]), Math.abs(expected) * 1e-13, "erfc(" + v[0] + ")");
		}
	}
}
