#!/usr/bin/env bash
# Builds hsdis, the plugin the JVM loads for -XX:+PrintAssembly and -XX:CompileCommand=print, so that compiled methods
# print as instructions instead of hex halfwords. bin/machcode.py reads either form -- without hsdis it disassembles
# the bytes with objdump -- so this is not required to annotate a profile; what it buys is reading the JVM's own output
# directly, which is how the architecture comparison workflow prints <system>-machcode.txt. Worth having before
# reasoning about what C2 did with a kernel: 20 calls in a hot method look alarming in hex and turn out, with
# mnemonics and the JVM's relocation comments beside them, to be uncommon traps.
# Usage: bin/build-hsdis.sh [--install] [--no-test]
#   --install  also copy the library into $JAVA_HOME/lib/server, where the JVM finds it without LD_LIBRARY_PATH.
#              Writes into the JDK installation; without it nothing outside target/ is touched.
#   --no-test  skip the self-test (a method compiled and printed, checked for instructions rather than hex).
#
# The backend is Capstone, which is one source file and one gcc command. The architecture comparison workflow builds
# the binutils backend instead, because it clones the JDK source there anyway; this is the same plugin interface and
# the one to use on a development machine, where cloning the JDK to disassemble a method is out of proportion.
#
# The sources are OpenJDK's own (src/utils/hsdis), downloaded from the tag of the JDK in JAVA_HOME so that the plugin
# interface matches the JVM that will load it, and kept in target/hsdis/src. Needs capstone's headers and library
# (Fedora: capstone-devel; Debian and derivatives: libcapstone-dev) and curl.
#
# Settings (environment): JAVA_HOME (the JDK to build for and to test with), OUT (target/hsdis).
set -eu

install=false
test=true
for arg in "$@"; do
	case "$arg" in
	--install) install=true ;;
	--no-test) test=false ;;
	*)
		echo "usage: bin/build-hsdis.sh [--install] [--no-test]" >&2
		exit 2
		;;
	esac
done

HERE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-$HERE/target/hsdis}"
: "${JAVA_HOME:?set JAVA_HOME to the JDK to build hsdis for}"

# The library's name is the one the JVM asks dlopen for: hsdis-<the JDK's name for this architecture>.so. The capstone
# architecture and mode are what the JDK's own build passes (make/Hsdis.gmk); capstone 6 renamed the AArch64 constant,
# so take whichever this header has.
case "$(uname -m)" in
x86_64) libarch=amd64 cs_arch=CS_ARCH_X86 cs_mode=CS_MODE_64 ;;
aarch64 | arm64) libarch=aarch64 cs_mode=CS_MODE_ARM ;;
*)
	echo "bin/build-hsdis.sh: no capstone architecture known for $(uname -m)" >&2
	exit 1
	;;
esac

capstone_flags="$(pkg-config --cflags --libs capstone 2>/dev/null || true)"
if [ -z "$capstone_flags" ]; then
	[ -f /usr/include/capstone/capstone.h ] ||
		{ echo "bin/build-hsdis.sh: capstone not found (Fedora: capstone-devel, Debian: libcapstone-dev)" >&2; exit 1; }
	capstone_flags="-I/usr/include/capstone -lcapstone"
fi
if [ "$libarch" = aarch64 ]; then
	header="$(pkg-config --variable=includedir capstone 2>/dev/null || echo /usr/include)/capstone/capstone.h"
	[ -f "$header" ] || header=/usr/include/capstone/capstone.h
	if grep -q 'CS_ARCH_AARCH64' "$header"; then cs_arch=CS_ARCH_AARCH64; else cs_arch=CS_ARCH_ARM64; fi
fi

# The JDK's own tag, as the workflow derives it: JAVA_RUNTIME_VERSION 27+35-2325 is the tag jdk-27+35. The update
# repository has the tags of the update releases and openjdk/jdk those of the feature releases, so try both.
runtime="$(sed -n 's/^JAVA_RUNTIME_VERSION="\(.*\)"$/\1/p' "$JAVA_HOME/release")"
feature="$(sed -n 's/^JAVA_VERSION="\([0-9]*\).*"$/\1/p' "$JAVA_HOME/release")"
tag="jdk-${runtime%%-*}"
tag_url="${tag//+/%2B}"

mkdir -p "$OUT/src"
if [ ! -f "$OUT/src/hsdis-capstone.c" ] || [ ! -f "$OUT/src/hsdis.h" ]; then
	got=
	for repo in "jdk${feature}u" jdk; do
		base="https://raw.githubusercontent.com/openjdk/$repo/$tag_url/src/utils/hsdis"
		if curl -fsS -o "$OUT/src/hsdis.h.part" "$base/hsdis.h" 2>/dev/null &&
			curl -fsS -o "$OUT/src/hsdis-capstone.c.part" "$base/capstone/hsdis-capstone.c" 2>/dev/null; then
			got="openjdk/$repo at $tag"
			break
		fi
	done
	if [ -z "$got" ]; then
		echo "bin/build-hsdis.sh: no hsdis sources for $tag; falling back to openjdk/jdk master" >&2
		base="https://raw.githubusercontent.com/openjdk/jdk/master/src/utils/hsdis"
		curl -fsS -o "$OUT/src/hsdis.h.part" "$base/hsdis.h"
		curl -fsS -o "$OUT/src/hsdis-capstone.c.part" "$base/capstone/hsdis-capstone.c"
		got="openjdk/jdk master"
	fi
	mv "$OUT/src/hsdis.h.part" "$OUT/src/hsdis.h"
	mv "$OUT/src/hsdis-capstone.c.part" "$OUT/src/hsdis-capstone.c"
	echo "sources: $got"
fi

lib="$OUT/hsdis-$libarch.so"
# shellcheck disable=SC2086
gcc -O2 -fPIC -shared -std=gnu11 \
	-D"CAPSTONE_ARCH=$cs_arch" -D"CAPSTONE_MODE=$cs_mode" \
	-I"$JAVA_HOME/include" -I"$JAVA_HOME/include/linux" \
	-o "$lib" "$OUT/src/hsdis-capstone.c" $capstone_flags
echo "built $lib ($cs_arch, $cs_mode, JDK $runtime)"

if $test; then
	# A method compiled and printed: with hsdis the JVM prints mnemonics, without it lines of hex halfwords.
	cat >"$OUT/HsdisSelfTest.java" <<'EOF'
public class HsdisSelfTest {
	static int work(int[] a) {
		int s = 0;
		for (int v : a) {
			s += v * 31;
		}
		return s;
	}

	public static void main(String[] args) {
		final int[] a = new int[1024];
		int s = 0;
		for (int r = 0; r < 20000; r++) {
			s += work(a);
		}
		if (s != 0) {
			throw new AssertionError();
		}
	}
}
EOF
	"$JAVA_HOME/bin/javac" -d "$OUT" "$OUT/HsdisSelfTest.java"
	LD_LIBRARY_PATH="$OUT${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$JAVA_HOME/bin/java" -Xbatch \
		-XX:+UnlockDiagnosticVMOptions -XX:CompileCommand=quiet -XX:CompileCommand=print,HsdisSelfTest::work \
		-cp "$OUT" HsdisSelfTest >"$OUT/selftest.log" 2>&1
	# The tab is the test: hsdis prints "mnemonic<tab><tab>operands", while the bytes the JVM prints without it are
	# space-separated halfwords, some of which ("c5f8", a VEX prefix) would pass for a mnemonic on their own.
	if grep -qE "^[[:space:]]+0x[0-9a-f]+:[[:space:]]+[a-z][a-z0-9.]*$(printf '\t')" "$OUT/selftest.log"; then
		echo "self-test: the JVM printed instructions, so it loaded $lib"
	else
		echo "self-test FAILED: the JVM printed no instructions; see $OUT/selftest.log" >&2
		exit 1
	fi
fi

if $install; then
	cp "$lib" "$JAVA_HOME/lib/server/"
	echo "installed $JAVA_HOME/lib/server/hsdis-$libarch.so; -XX:+PrintAssembly now needs no LD_LIBRARY_PATH"
else
	echo "use it with:  LD_LIBRARY_PATH=$OUT <java> -XX:+UnlockDiagnosticVMOptions -XX:PrintAssemblyOptions=intel \\"
	echo "                -XX:CompileCommand=quiet -XX:CompileCommand=print,pkg.Class::method ..."
	echo "or run with --install to put it in the JDK, where the JVM finds it by itself."
fi
