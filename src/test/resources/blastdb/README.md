# Test BLAST databases

Made with BLAST+ 2.17.0. Only the `.pin`, `.psq` and `.phr` files (and the `.pal` alias file) are kept; blastp4j
does not read the others.

    makeblastdb -in ids.fa -dbtype prot -parse_seqids -out ids_p
    makeblastdb -in ids.fa -dbtype prot -out ids_n
    zcat ../seg/sprot_sample.fasta.gz > sample.fa
    makeblastdb -in sample.fa -dbtype prot -parse_seqids -out multi -max_file_sz 30KB

`ids.fa` has one sequence for each kind of identifier; `q.fa` is the query. The expected outputs:

    blastp -query q.fa -db ids_p -evalue 1000 -outfmt "6 qseqid sseqid sacc saccver stitle slen evalue bitscore" > ids_p.tsv
    blastp -query q.fa -db ids_n -evalue 1000 -outfmt "6 qseqid sseqid sacc saccver stitle slen evalue bitscore" > ids_n.tsv

`multi` is the 400-sequence sample split into six volumes; blastp against it gives `../cli/sample_std_btop.tsv.gz`.
