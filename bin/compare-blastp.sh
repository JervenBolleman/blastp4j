#!/usr/bin/env bash
# Compares BLAST+ blastp, blastp4j on the JVM and a native blastp4j on one search.
# Usage: bin/compare-blastp.sh [all|blast|jvm|native]...   (default: all)
#
# 1. ROUNDS timing rounds (default 5), the chosen systems interleaved in each round (/usr/bin/time: wall, user, sys).
# 2. One perf stat run per system (instructions, cycles, branches, branch misses), after all timing rounds.
# 3. For the JVM, one run with async-profiler (CPU samples, flame graph in jvm-profile.html).
# Every output is compared with the first one written (BLAST+'s when blast is chosen).
#
# Settings (environment): QUERY, DB, THREADS (8), ROUNDS (5), OUT (target/compare), BLASTP (blastp on PATH),
# JAVA_HOME or JAVA, NATIVE (the native executable, needed for native), ASYNC_PROFILER (libasyncProfiler.so; without
# it the profile run is skipped). The JVM runs the jar built by "mvn package".
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

SYSTEMS=()
for a in "${@:-all}"; do
	case "$a" in
	all) SYSTEMS=(blast jvm native) ;;
	blast | jvm | native) SYSTEMS+=("$a") ;;
	*)
		echo "usage: $0 [all|blast|jvm|native]..." >&2
		exit 1
		;;
	esac
done

JAR="$(ls "$HERE"/target/blastp4j-*.jar 2>/dev/null | grep -v -e sources -e javadoc | head -1)"
for s in "${SYSTEMS[@]}"; do
	case $s in
	jvm) [ -n "$JAR" ] || { echo "compare-blastp: no jar in $HERE/target; run 'mvn package' first" >&2; exit 1; } ;;
	native) [ -x "${NATIVE:-}" ] || { echo "compare-blastp: set NATIVE to the native blastp4j executable" >&2; exit 1; } ;;
	esac
done
mkdir -p "$OUT"

# The command of a system, before the blastp options; extra JVM options as further arguments.
command_of() {
	case $1 in
	blast) echo "$BLASTP" ;;
	jvm) echo "$JAVA ${2:-} -jar $JAR" ;;
	native) echo "$NATIVE" ;;
	esac
}

REFERENCE=""
# Runs one search, its output in $OUT/$2.tsv, and compares it with the reference.
search() {
	local system=$1 name=$2
	shift 2
	"$@" $(command_of "$system" "${JVM_EXTRA:-}") -query "$QUERY" -db "$DB" -outfmt '6 std btop' -num_threads "$THREADS" \
		-out "$OUT/$name.tsv" >"$OUT/$name.log" 2>&1
	if [ -z "$REFERENCE" ]; then
		REFERENCE="$OUT/$name.tsv"
		SAME="reference"
	elif cmp -s "$REFERENCE" "$OUT/$name.tsv"; then
		SAME="identical"
	else
		SAME="DIFFERENT from $(basename "$REFERENCE")"
	fi
}

echo "systems: ${SYSTEMS[*]}; query $QUERY; db $DB; $THREADS threads; output in $OUT"
echo "== timing, $ROUNDS rounds (wall s, user s, sys s)"
for r in $(seq 1 "$ROUNDS"); do
	for s in "${SYSTEMS[@]}"; do
		search "$s" "time-$s-$r" /usr/bin/time -o "$OUT/time-$s-$r.txt" -f "%e %U %S"
		echo "round $r $s: $(cat "$OUT/time-$s-$r.txt") ($SAME)"
	done
done
for s in "${SYSTEMS[@]}"; do
	echo "$s median wall: $(cat "$OUT"/time-"$s"-*.txt | cut -d' ' -f1 | sort -n | awk '{a[NR]=$1} END{print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2) " s"}')"
done

echo "== perf stat"
for s in "${SYSTEMS[@]}"; do
	search "$s" "perf-$s" perf stat -x, -e instructions:u,cycles:u,branches:u,branch-misses:u -o "$OUT/perf-$s.txt" --
	echo "$s: $(awk -F, 'NR>1 && $1 ~ /^[0-9]/{printf "%s=%.1fG ", $3, $1/1e9}' "$OUT/perf-$s.txt")($SAME)"
done

if [[ " ${SYSTEMS[*]} " == *" jvm "* ]]; then
	echo "== async-profiler (JVM)"
	if [ -f "$ASYNC_PROFILER" ]; then
		JVM_EXTRA="-agentpath:$ASYNC_PROFILER=start,event=cpu,file=$OUT/jvm-profile.html" search jvm "profile-jvm"
		echo "flame graph: $OUT/jvm-profile.html ($SAME)"
	else
		echo "skipped: set ASYNC_PROFILER to libasyncProfiler.so"
	fi
fi
