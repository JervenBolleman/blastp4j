#!/usr/bin/env bash
# Regenerates the BLAST+ reference results in src/test/resources/reference.
#
#   BLAST_BIN=~/opt/ncbi-blast-2.17.0+/bin SPROT=~/data/uniprot reference/make_reference.sh
#
# Needs BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and UniProtKB/Swiss-Prot
# (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz).
# The database is built with makeblastdb if it does not exist. Single-threaded, so the times are comparable.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BLAST_BIN="${BLAST_BIN:-$HOME/opt/ncbi-blast-2.17.0+/bin}"
SPROT="${SPROT:-$HOME/data/uniprot}"
OUT="$ROOT/src/test/resources/reference"
QUERIES="$ROOT/src/test/resources/queries/uniprot_human.fasta"
DB="$SPROT/blastdb/sprot"
FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score qlen slen btop"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [[ ! -f "$DB.pin" ]]; then
	[[ -f "$SPROT/uniprot_sprot.fasta" ]] || gunzip -k "$SPROT/uniprot_sprot.fasta.gz"
	mkdir -p "$SPROT/blastdb"
	"$BLAST_BIN/makeblastdb" -in "$SPROT/uniprot_sprot.fasta" -dbtype prot -parse_seqids -out "$DB" \
		-title "UniProtKB/Swiss-Prot"
fi

"$BLAST_BIN/blastp" -version
for cbs in 2 0; do
	/usr/bin/time -f "comp_based_stats $cbs: %e s wall, %M KB max RSS" \
		"$BLAST_BIN/blastp" -query "$QUERIES" -db "$DB" -num_threads 1 -comp_based_stats "$cbs" \
		-outfmt "$FORMAT" -out "$WORK/out.tsv"
	gzip -9nc "$WORK/out.tsv" > "$OUT/human_vs_sprot.cbs$cbs.tsv.gz"
	echo "human_vs_sprot.cbs$cbs.tsv.gz: $(wc -l < "$WORK/out.tsv") hits"
done

# Ungapped search: its alignments are the ungapped extensions of the seeding stage.
/usr/bin/time -f "ungapped: %e s wall, %M KB max RSS" \
	"$BLAST_BIN/blastp" -query "$QUERIES" -db "$DB" -num_threads 1 -comp_based_stats 0 -ungapped \
	-outfmt "$FORMAT" -out "$WORK/out.tsv"
gzip -9nc "$WORK/out.tsv" > "$OUT/human_vs_sprot.ungapped.tsv.gz"
echo "human_vs_sprot.ungapped.tsv.gz: $(wc -l < "$WORK/out.tsv") hits"
