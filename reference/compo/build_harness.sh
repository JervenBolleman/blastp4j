#!/usr/bin/env bash
# Builds harness.c against BLAST+'s composition adjustment library and runs it on a cases file:
#
#   reference/compo/build_harness.sh cases.txt > blast.txt
#
# The library sources (public domain, NCBI C++ Toolkit) are downloaded into a temporary directory. The harness reads
# the cases written by blastp4j (compositions as hex doubles) and prints, for each, the status, matrix adjustment rule
# and lambda ratio of Blast_AdjustScores and the adjusted 28 x 28 matrix; src/test/resources/compo holds such a pair.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CASES="$(realpath "$1")"
REPO=https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
mkdir -p inc/algo/blast/composition_adjustment inc/algo/blast/core
for f in composition_adjustment optimize_target_freq nlm_linear_algebra matrix_frequency_data compo_mode_condition \
	unified_pvalues; do
	curl -sfL -o $f.c $REPO/src/algo/blast/composition_adjustment/$f.c
done
for f in composition_constants composition_adjustment compo_heap compo_mode_condition matrix_frequency_data \
	nlm_linear_algebra optimize_target_freq redo_alignment smith_waterman unified_pvalues; do
	curl -sfL -o inc/algo/blast/composition_adjustment/$f.h $REPO/include/algo/blast/composition_adjustment/$f.h
done
for f in ncbi_std blast_export ncbi_math blast_def; do
	curl -sfL -o inc/algo/blast/core/$f.h $REPO/include/algo/blast/core/$f.h
done
curl -sfL -o matrix_freq_ratios.c $REPO/src/algo/blast/core/matrix_freq_ratios.c
curl -sfL -o blast_stat.c $REPO/src/algo/blast/core/blast_stat.c
# The toolkit's basic types, without the rest of the toolkit.
cat > inc/algo/blast/core/blast_toolkit.h <<'H'
#ifndef _BLAST_TOOLKIT__H
#define _BLAST_TOOLKIT__H
#include <stdint.h>
typedef int8_t Int1; typedef uint8_t Uint1; typedef int16_t Int2; typedef uint16_t Uint2;
typedef int32_t Int4; typedef uint32_t Uint4; typedef int64_t Int8; typedef uint64_t Uint8;
#endif
H
# BLOSUM62 frequency ratios and BLAST's Newton iteration for lambda, copied from the sources.
python3 -I - matrix_freq_ratios.c > freqratios.inc <<'P'
import sys
s = open(sys.argv[1]).read()
i = s.index('{', s.index('BLOSUM62_FREQRATIOS[BLASTAA_SIZE][BLASTAA_SIZE] ='))
d = 0
for j in range(i, len(s)):
    d += {'{': 1, '}': -1}.get(s[j], 0)
    if d == 0:
        break
print('static const double FREQRATIOS[28][28] = ' + s[i:j + 1] + ';')
P
sed -n "$(grep -n '^NlmKarlinLambdaNR' blast_stat.c | cut -d: -f1),+73p" blast_stat.c > nr.inc
cp "$HERE/harness.c" .
gcc -O2 -Iinc -o harness harness.c composition_adjustment.c optimize_target_freq.c nlm_linear_algebra.c \
	matrix_frequency_data.c compo_mode_condition.c unified_pvalues.c -lm
./harness < "$CASES"
