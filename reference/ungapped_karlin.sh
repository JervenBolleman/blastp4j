#!/usr/bin/env bash
# Records blastp's ungapped Karlin-Altschul parameters (the first "Lambda K H" line of the text report) for each
# reference query, into src/test/resources/reference/ungapped_karlin.tsv. They depend only on the query composition,
# so each query is searched against itself.
#
#   BLAST_BIN=~/opt/ncbi-blast-2.17.0+/bin reference/ungapped_karlin.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BLAST_BIN="${BLAST_BIN:-$HOME/opt/ncbi-blast-2.17.0+/bin}"
QUERIES="$ROOT/src/test/resources/queries/uniprot_human.fasta"
OUT="$ROOT/src/test/resources/reference/ungapped_karlin.tsv"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

awk -v dir="$WORK" '/^>/ { n++; f = sprintf("%s/q%02d.fa", dir, n) } { print > f }' "$QUERIES"
printf 'qseqid\tlambda\tK\tH\n' > "$OUT"
for q in "$WORK"/q*.fa; do
	id=$(head -1 "$q" | cut -c2- | cut -d' ' -f1)
	"$BLAST_BIN/blastp" -query "$q" -subject "$q" -comp_based_stats 0 -outfmt 0 |
		awk -v id="$id" '/^Lambda/ && !done { getline; printf "%s\t%s\t%s\t%s\n", id, $1, $2, $3; done = 1 }' >> "$OUT"
done
echo "$OUT: $(($(wc -l < "$OUT") - 1)) queries"
