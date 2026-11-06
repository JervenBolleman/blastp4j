package swiss.sib.blastp4j.util;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import jdk.incubator.vector.IntVector;

/**
 * Whether the kernels that are otherwise at most 256 bits wide use the machine's 512-bit vectors: on AMD Zen 5 (CPU
 * family 26) with AVX-512, whose vector units are 512 bits wide (Zen 4 splits each 512-bit operation in two), and only
 * where the JVM's preferred vectors are 512 bits. {@code -Dblastp4j.wideVectors=true} or {@code false} overrides the
 * detection. A constant, read once.
 */
public final class WideVectors {
	public static final boolean ENABLED = enabled();

	private WideVectors() {
	}

	private static boolean enabled() {
		final String forced = System.getProperty("blastp4j.wideVectors");
		if (forced != null) {
			return Boolean.parseBoolean(forced);
		}
		return IntVector.SPECIES_PREFERRED.length() >= 16 && zen5();
	}

	/** Whether /proc/cpuinfo names an AMD CPU of family 26 with AVX-512 (Linux only; elsewhere false). */
	static boolean zen5() {
		try {
			return zen5(Files.readAllLines(Path.of("/proc/cpuinfo")));
		} catch (IOException | RuntimeException e) {
			return false;
		}
	}

	/** {@link #zen5()} for the lines of a /proc/cpuinfo: those of its first CPU. */
	static boolean zen5(List<String> cpuinfo) {
		boolean amd = false, family26 = false, avx512 = false;
		for (String line : cpuinfo) {
			if (line.isBlank()) {
				break;
			}
			final int colon = line.indexOf(':');
			if (colon < 0) {
				continue;
			}
			final String key = line.substring(0, colon).trim(), value = line.substring(colon + 1).trim();
			switch (key) {
			case "vendor_id" -> amd = value.equals("AuthenticAMD");
			case "cpu family" -> family26 = value.equals("26");
			case "flags" -> avx512 = (" " + value + " ").contains(" avx512f ");
			default -> {
			}
			}
		}
		return amd && family26 && avx512;
	}
}
