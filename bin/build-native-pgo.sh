#!/usr/bin/env sh
# Builds a native blastp4j with profile-guided optimisation (GraalVM native-image, -O3): an instrumented build, a training
# run that writes default.iprof, then the optimised build. Needs "mvn compile" to have run (target/classes and target/cp.txt).
# Usage: bin/build-native-pgo.sh <training.fasta> <blast db prefix> <output directory>
# NATIVE_IMAGE names the GraalVM native-image to use (default: native-image on PATH). The Vector API module and the
# build-time initialization of the word finder come from META-INF/native-image in target/classes. The default launcher
# (bin/blastp4j) is unchanged.
set -e
if [ $# -ne 3 ]; then
	echo "usage: $0 <training.fasta> <blast db prefix> <output directory>" >&2
	exit 1
fi
HERE="$(cd "$(dirname "$0")/.." && pwd)"
TRAIN="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
DB="$2"
OUT="$(mkdir -p "$3" && cd "$3" && pwd)"
CLASSES="$HERE/target/classes"
CP_FILE="$HERE/target/cp.txt"
if [ ! -d "$CLASSES" ] || [ ! -f "$CP_FILE" ]; then
	echo "build-native-pgo: run 'mvn compile' first (needs $CLASSES and $CP_FILE)" >&2
	exit 1
fi
NI="${NATIVE_IMAGE:-native-image}"
CP="$CLASSES:$(cat "$CP_FILE")"
OPTS="-O3"

echo "1/3 instrumented build"
"$NI" $OPTS --pgo-instrument -cp "$CP" swiss.sib.blastp4j.cli.Blastp -o "$OUT/blastp4j-instrumented"

echo "2/3 training run (writes default.iprof in the working directory)"
(cd "$OUT" && ./blastp4j-instrumented -query "$TRAIN" -db "$DB" -outfmt '6 std btop' -num_threads 8 -out "$OUT/training.tsv")

echo "3/3 optimised build with the profile"
"$NI" $OPTS --pgo="$OUT/default.iprof" -cp "$CP" swiss.sib.blastp4j.cli.Blastp -o "$OUT/blastp4j-pgo"
echo "built $OUT/blastp4j-pgo"
