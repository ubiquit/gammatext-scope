#!/usr/bin/env bash
# GENCODE GFF3 -> BED6 files for paftools liftover, in three groups, each with a key map.
#
# paftools liftover names every lifted interval by its source coordinates
# (chr_start_end), so the attributes of a feature have to be looked up by those
# coordinates afterwards, and features that share a span collide. A gene whose
# span equals one of its transcripts would lose its "gene" line. The features are
# therefore lifted in three groups: protein-coding genes, other genes and
# everything else, so a gene can only collide with another gene of the same
# group and the same span (two protein-coding genes with identical spans in
# GENCODE v44, which does not change the gene mask).
#
# Chromosome names are changed from chr1..chrY, chrM to the Ensembl-style
# 1..Y, MT used by the reference FASTA.
#
# usage: inputs/gencode_to_groups.sh gencode.v44.annotation.gff3.gz <outdir>
# writes <outdir>/{pcgene,othergene,rest}.bed6 and .map.tsv (<chr_start_end>\t<type>\t<source>\t<attributes>)
set -euo pipefail
GFF=${1:?gencode gff3(.gz)}; OUT=${2:?outdir}
mkdir -p "$OUT"
{ case "$GFF" in *.gz) zcat "$GFF";; *) cat "$GFF";; esac; } | awk -F'\t' -v OFS='\t' -v dir="$OUT" '
  /^#/ { next }
  NF >= 9 {
    chr = $1; sub(/^chr/, "", chr); if (chr == "M") chr = "MT";
    g = ($3 == "gene") ? (($9 ~ /gene_type=protein_coding/) ? "pcgene" : "othergene") : "rest";
    print chr, $4 - 1, $5, chr "_" ($4 - 1) "_" $5, 0, $7 > (dir "/" g ".bed6");
    print chr "_" ($4 - 1) "_" $5, $3, $2, $9 > (dir "/" g ".map.tsv");
  }'
for g in pcgene othergene rest; do echo "$g: $(wc -l < "$OUT/$g.bed6") features"; done
