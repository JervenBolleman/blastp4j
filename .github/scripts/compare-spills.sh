#!/usr/bin/env bash
# Compares register spills in the X-drop alignment loop: BLAST+'s Blast_SemiGappedAlign (gcc) against blastp4j's
# GappedAligner.semiGapped (C2). Runs one search with each under "perf record -e cpu-clock" (a software clock, so it
# also works in virtual machines without hardware counters), prints C2's assembly with hsdis, and counts the stack
# accesses in each function: statically, and weighted by the samples. Prints a Markdown table.
#
# Settings (environment): QUERY, DB, THREADS (nproc), OUT (target/spills), BLASTP (blastp on PATH), JAVA_HOME,
# BLASTP4J_OPTS (extra JVM options). hsdis must be where the JVM finds it (lib/server of the JDK, or LD_LIBRARY_PATH).
set -e
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
QUERY="${QUERY:-$HERE/src/test/resources/queries/uniprot_human.fasta}"
DB="${DB:?set DB to a BLAST protein database}"
THREADS="${THREADS:-$(nproc)}"
OUT="${OUT:-$HERE/target/spills}"
BLASTP="$(command -v "${BLASTP:-blastp}")"
COUNT="$HERE/.github/scripts/loop-stack-accesses.py"
METHOD=swiss.sib.blastp4j.align.GappedAligner::semiGapped
mkdir -p "$OUT"

perf record -q -e cpu-clock -F 4999 -o "$OUT/blast.perf" -- \
	"$BLASTP" -query "$QUERY" -db "$DB" -outfmt '6 std btop' -num_threads "$THREADS" -out "$OUT/blast.tsv"
perf script -F ip -i "$OUT/blast.perf" >"$OUT/blast.ips" 2>/dev/null
objdump -d --no-show-raw-insn --disassemble=Blast_SemiGappedAlign "$BLASTP" >"$OUT/Blast_SemiGappedAlign.s"

BLASTP4J_OPTS="$BLASTP4J_OPTS -XX:+UnlockDiagnosticVMOptions -XX:CompileCommand=quiet -XX:CompileCommand=print,$METHOD" \
	perf record -q -e cpu-clock -F 4999 -o "$OUT/jvm.perf" -- \
	"$HERE/bin/blastp4j" -query "$QUERY" -db "$DB" -outfmt '6 std btop' -num_threads "$THREADS" \
	-out "$OUT/jvm.tsv" >"$OUT/semiGapped.s" 2>"$OUT/jvm.log"
perf script -F ip -i "$OUT/jvm.perf" >"$OUT/jvm.ips" 2>/dev/null
grep -q '^\[Verified Entry Point\]' "$OUT/semiGapped.s" ||
	{ echo "compare-spills: no disassembly in $OUT/semiGapped.s; is hsdis installed?" >&2; exit 1; }

if cmp -s "$OUT/blast.tsv" "$OUT/jvm.tsv"; then SAME="identical"; else SAME="DIFFERENT"; fi
echo "Search outputs: $SAME. Stack accesses in the X-drop alignment function ($(uname -m)):"
echo
python3 "$COUNT" --header
python3 "$COUNT" --objdump "$OUT/Blast_SemiGappedAlign.s" --samples "$OUT/blast.ips" \
	--label "BLAST+ Blast_SemiGappedAlign ($(readelf -p .comment "$BLASTP" 2>/dev/null | grep -o 'GCC: ([^)]*) [0-9.]*' | head -1))"
python3 "$COUNT" --hotspot "$OUT/semiGapped.s" --method semiGapped --samples "$OUT/jvm.ips" \
	--label "blastp4j semiGapped (C2, $("${JAVA_HOME:+$JAVA_HOME/bin/}java" -version 2>&1 | head -1 | cut -d'"' -f2))"
