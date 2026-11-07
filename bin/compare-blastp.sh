#!/usr/bin/env bash
# Compares BLAST+ blastp, blastp4j on the JVM and a native blastp4j on one search.
# Usage: bin/compare-blastp.sh [all|blast|jvm|ref|almost|native|blastfast|fast]...   (default: all, which is blast,
# jvm and native)
# almost is almostblastp (cli.AlmostBlastp) from the same jar as jvm, whose results may differ slightly from blastp's.
# blastfast is BLAST+ blastp -task blastp-fast, and fast fastblastp (cli.FastBlastp) from the same jar, which gives
# blastp-fast's results.
#
# 1. ROUNDS timing rounds (default 5), the chosen systems interleaved in each round (/usr/bin/time: wall, user, sys).
# 2. One perf stat run per system (instructions, cycles, branches, branch misses), after all timing rounds; skipped
#    where perf has no hardware counters.
# 3. For the JVM, one run with async-profiler (CPU samples, flame graph in jvm-profile.html), and one that prints the
#    machine code of the seed and align packages with the samples of each instruction (machcode.txt).
# Every output is compared with the first one written of its task (BLAST+'s when blast, or blastfast, is chosen):
# blastfast and fast with each other, the others with each other.
#
# Settings (environment): QUERY, DB, THREADS (8), ROUNDS (5), OUT (target/compare), BLASTP (blastp on PATH),
# TIME (GNU time, for -o and -f; default /usr/bin/time, on macOS gtime from Homebrew's gnu-time), JAVA_HOME or JAVA, NATIVE (the native executable, needed for native), ASYNC_PROFILER (libasyncProfiler.so; without
# it the profile run is skipped). The JVM runs the jar built by "mvn package", with this machine's options from
# bin/jvm.options (as bin/blastp4j). ref runs another build's jar, REF_JAR, with the options of that build's tree,
# REF_HOME (default this tree's): blastp4j at another commit, to compare code and options on the same machine.
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
QUERY="${QUERY:-$HERE/src/test/resources/queries/uniprot_human.fasta}"
DB="${DB:-$HOME/data/uniprot/blastdb/sprot}"
THREADS="${THREADS:-8}"
ROUNDS="${ROUNDS:-5}"
OUT="${OUT:-$HERE/target/compare}"
BLASTP="${BLASTP:-blastp}"
JAVA="${JAVA:-${JAVA_HOME:+$JAVA_HOME/bin/}java}"
ASYNC_PROFILER="${ASYNC_PROFILER:-}"
TIME="${TIME:-/usr/bin/time}"

SYSTEMS=()
for a in "${@:-all}"; do
	case "$a" in
	all) SYSTEMS=(blast jvm native) ;;
	blast | jvm | ref | almost | native | blastfast | fast) SYSTEMS+=("$a") ;;
	*)
		echo "usage: $0 [all|blast|jvm|ref|almost|native|blastfast|fast]..." >&2
		exit 1
		;;
	esac
done

JAR="$(ls "$HERE"/target/blastp4j-*.jar 2>/dev/null | grep -v -e sources -e javadoc | head -1)"
for s in "${SYSTEMS[@]}"; do
	case $s in
	jvm | almost | fast) [ -n "$JAR" ] || { echo "compare-blastp: no jar in $HERE/target; run 'mvn package' first" >&2; exit 1; } ;;
	ref) [ -f "${REF_JAR:-}" ] || { echo "compare-blastp: set REF_JAR to the jar to compare with" >&2; exit 1; } ;;
	native) [ -x "${NATIVE:-}" ] || { echo "compare-blastp: set NATIVE to the native blastp4j executable" >&2; exit 1; } ;;
	esac
done
mkdir -p "$OUT"

# The command of a system, before the blastp options; extra JVM options as further arguments.
command_of() {
	case $1 in
	blast) echo "$BLASTP" ;;
	jvm) echo "$JAVA $("$HERE/bin/jvm-options") ${2:-} -jar $JAR" ;;
	almost) echo "$JAVA $("$HERE/bin/jvm-options") ${2:-} -cp $JAR swiss.sib.blastp4j.cli.AlmostBlastp" ;;
	blastfast) echo "$BLASTP -task blastp-fast" ;;
	fast) echo "$JAVA $("$HERE/bin/jvm-options") ${2:-} -cp $JAR swiss.sib.blastp4j.cli.FastBlastp" ;;
	ref) echo "$JAVA $("${REF_HOME:-$HERE}/bin/jvm-options") ${2:-} -jar $REF_JAR" ;;
	native) echo "$NATIVE" ;;
	esac
}

# The reference output of each task: blastp's, and blastp -task blastp-fast's (plain variables: macOS has bash 3).
REFERENCE_blastp=""
REFERENCE_fast=""
# Runs one search, its output in $OUT/$2.tsv, and compares it with its task's reference.
search() {
	local system=$1 name=$2 task=blastp reference
	case $1 in blastfast | fast) task=fast ;; esac
	shift 2
	"$@" $(command_of "$system" "${JVM_EXTRA:-}") -query "$QUERY" -db "$DB" -outfmt '6 std btop' -num_threads "$THREADS" \
		-out "$OUT/$name.tsv" >"$OUT/$name.log" 2>&1
	eval "reference=\$REFERENCE_$task"
	if [ -z "$reference" ]; then
		eval "REFERENCE_$task=\$OUT/\$name.tsv"
		SAME="reference"
	elif cmp -s "$reference" "$OUT/$name.tsv"; then
		SAME="identical"
	else
		SAME="DIFFERENT from $(basename "$reference")"
	fi
}

echo "systems: ${SYSTEMS[*]}; query $QUERY; db $DB; $THREADS threads; output in $OUT"
echo "== timing, $ROUNDS rounds (wall s, user s, sys s)"
for r in $(seq 1 "$ROUNDS"); do
	for s in "${SYSTEMS[@]}"; do
		search "$s" "time-$s-$r" "$TIME" -o "$OUT/time-$s-$r.txt" -f "%e %U %S"
		echo "round $r $s: $(cat "$OUT/time-$s-$r.txt") ($SAME)"
	done
done
for s in "${SYSTEMS[@]}"; do
	echo "$s median wall: $(cat "$OUT"/time-"$s"-*.txt | cut -d' ' -f1 | sort -n | awk '{a[NR]=$1} END{print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2) " s"}')"
done

echo "== perf stat"
# Virtual machines often have no hardware counters; perf then reports them as not supported.
if perf stat -x, -e cycles:u -o "$OUT/perf-check.txt" -- true 2>/dev/null && grep -q '^[0-9]' "$OUT/perf-check.txt"; then
	for s in "${SYSTEMS[@]}"; do
		search "$s" "perf-$s" perf stat -x, -e instructions:u,cycles:u,branches:u,branch-misses:u -o "$OUT/perf-$s.txt" --
		echo "$s: $(awk -F, 'NR>1 && $1 ~ /^[0-9]/{printf "%s=%.1fG ", $3, $1/1e9}' "$OUT/perf-$s.txt")($SAME)"
	done
else
	echo "skipped: no hardware counters (perf stat -e cycles:u)"
fi

if [[ " ${SYSTEMS[*]} " == *" jvm "* ]]; then
	echo "== async-profiler (JVM)"
	if [ -f "$ASYNC_PROFILER" ]; then
		JVM_EXTRA="-agentpath:$ASYNC_PROFILER=start,event=cpu,file=$OUT/jvm-profile.jfr" search jvm "profile-jvm"
		# jfrconv is in the bin directory next to the library's lib directory.
		JFRCONV="$(dirname "$ASYNC_PROFILER")/../bin/jfrconv"
		"$JFRCONV" --norm -o html "$OUT/jvm-profile.jfr" "$OUT/jvm-profile.html" >/dev/null
		echo "flame graph: $OUT/jvm-profile.html ($SAME)"
		# The methods with the most samples of their own (the top frame), hidden-class copies merged by --norm.
		"$JFRCONV" --norm --simple -o collapsed "$OUT/jvm-profile.jfr" "$OUT/jvm-profile.collapsed" >/dev/null
		echo "self time, top 15 methods:"
		awk '{ n = $NF; sub(/ [0-9]+$/, ""); k = split($0, f, ";"); self[f[k]] += n; total += n }
			END { for (m in self) printf "%6.2f%%  %s\n", 100 * self[m] / total, m }' "$OUT/jvm-profile.collapsed" |
			sort -rn | head -15
	else
		echo "skipped: set ASYNC_PROFILER to libasyncProfiler.so"
	fi

	echo "== machine code (JVM)"
	# The compiled code of the seeding and alignment packages, with the samples of the same run per instruction
	# (bin/machcode.py): cycles where perf has hardware counters, a software clock otherwise. With hsdis the JVM prints
	# the disassembly, without it the bytes, which objdump disassembles.
	if command -v perf >/dev/null && command -v objdump >/dev/null; then
		EVENT=cpu-clock:u
		grep -q '^[0-9]' "$OUT/perf-check.txt" 2>/dev/null && EVENT=cycles:u
		PRINT="-XX:+UnlockDiagnosticVMOptions -XX:CompileCommand=quiet"
		PRINT+=" -XX:CompileCommand=print,swiss.sib.blastp4j.seed.*::* -XX:CompileCommand=print,swiss.sib.blastp4j.align.*::*"
		[ "$(uname -m)" = x86_64 ] && PRINT+=" -XX:PrintAssemblyOptions=intel"
		JVM_EXTRA="$PRINT" search jvm "machcode-jvm" perf record -q -e "$EVENT" -F 4999 -o "$OUT/machcode.perf" --
		perf script -F ip -i "$OUT/machcode.perf" >"$OUT/machcode.ips" 2>/dev/null
		python3 "$HERE/bin/machcode.py" "$OUT/machcode-jvm.log" "$OUT/machcode.ips" >"$OUT/machcode.txt"
		echo "$EVENT samples ($SAME); the full listing per method in machcode.txt:"
		sed -n '/^==/q;p' "$OUT/machcode.txt"
	else
		echo "skipped: needs perf and objdump"
	fi
fi
