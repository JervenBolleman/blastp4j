package swiss.sib.blastp4j.util;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

class WideVectorsTest {
	private static List<String> cpu(String vendor, int family, String flags) {
		return List.of("processor\t: 0", "vendor_id\t: " + vendor, "cpu family\t: " + family,
				"model name\t: some CPU", "flags\t\t: " + flags, "", "processor\t: 1", "vendor_id\t: GenuineIntel");
	}

	@Test
	void zen5WithAvx512() {
		assertTrue(WideVectors.zen5(cpu("AuthenticAMD", 26, "fpu sse2 avx2 avx512f avx512bw")));
	}

	@Test
	void notZen4NorIntelNorWithoutAvx512() {
		assertFalse(WideVectors.zen5(cpu("AuthenticAMD", 25, "fpu avx2 avx512f")));
		assertFalse(WideVectors.zen5(cpu("GenuineIntel", 26, "fpu avx2 avx512f")));
		assertFalse(WideVectors.zen5(cpu("AuthenticAMD", 26, "fpu avx2 avx512fx")));
		assertFalse(WideVectors.zen5(List.of()));
	}
}
