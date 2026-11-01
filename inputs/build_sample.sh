#!/usr/bin/env bash
# Build one sample of the data layout from a reference genome and a VCF:
#   <data-root>/<Sample>/chr<N>.<id>.hap1.txt           consensus sequence per chromosome, no FASTA header
#   <data-root>/<Sample>_AccessionProcessing/genomic.gff annotation lifted onto the consensus (GFF3)
#   <data-root>/<Sample>_AccessionProcessing/<Sample>_Chromosomes.json
# after which `gammatext-scope prepare-annotation --sample <Sample>` produces <Sample>_GeneInfo.
#
# usage: inputs/build_sample.sh <Sample> <id> <ref.fa> <vcf-or-pattern> "<bcftools consensus options>" <groups-dir> <paftools.js>
#   <Sample>          sample folder name, e.g. HomoSapiens_BreastCancer_Normal_v2
#   <id>              sequence id used in the file names, e.g. HCC1395BL
#   <ref.fa>          GRCh38 primary assembly, sequences named 1..22 X Y MT, with .fai
#   <vcf-or-pattern>  one bgzipped+indexed VCF for all chromosomes, or a per-chromosome
#                     pattern with {c} for the chromosome name, e.g. vcf/{c}.normal.vcf.gz
#   options           passed to bcftools consensus, e.g. "-s COLO829 -H A" (the ALT allele of
#                     every genotype of that sample) or "" for a sites-only VCF
#   <groups-dir>      the three BED6 files and key maps written by inputs/gencode_to_groups.sh
#   <paftools.js>     from minimap2 (run with k8)
# environment: GAMMATEXT_DATA (data root, default ./data), WORK (default ./work), THREADS (default 8)
# tools: bcftools, samtools, minimap2, k8, awk, sort
set -euo pipefail
SAMPLE=${1:?sample}; ID=${2:?id}; REF=${3:?ref.fa}; VCFPAT=${4:?vcf}; OPTS=${5-}; GROUPS=${6:?groups-dir}; PAFTOOLS=${7:?paftools.js}
DATA=${GAMMATEXT_DATA:-data}; WORK=${WORK:-work}/$SAMPLE; THREADS=${THREADS:-8}
CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
mkdir -p "$WORK/chr" "$DATA/$SAMPLE" "$DATA/${SAMPLE}_AccessionProcessing" "$DATA/${SAMPLE}_GeneInfo"

echo "=== 1. consensus per chromosome"
for c in $CHRS; do
  fa="$WORK/chr/chr$c.$ID.hap1.fa"
  [ -s "$fa" ] && continue
  vcf="${VCFPAT//\{c\}/$c}"
  # shellcheck disable=SC2086
  samtools faidx "$REF" "$c" | bcftools consensus $OPTS "$vcf" > "$fa" 2> "$fa.log"
  grep -q "does not match" "$fa.log" && echo "   WARNING chr$c: REF mismatches, see $fa.log" || true
  grep -v '^>' "$fa" > "$DATA/$SAMPLE/chr$c.$ID.hap1.txt"
done
files=(); for c in $CHRS; do files+=("$WORK/chr/chr$c.$ID.hap1.fa"); done
cat "${files[@]}" > "$WORK/$SAMPLE.consensus.fa"; samtools faidx "$WORK/$SAMPLE.consensus.fa"
echo "   $(ls "$DATA/$SAMPLE" | wc -l) chromosome files, $(awk '{s+=$2}END{print s}' "$WORK/$SAMPLE.consensus.fa.fai") bp"

echo "=== 2. alignment consensus <- reference (the reference is the query)"
if [ ! -s "$WORK/from_REF.sorted.paf" ]; then
  minimap2 -t "$THREADS" -x asm5 -c --secondary=no "$WORK/$SAMPLE.consensus.fa" "$REF" > "$WORK/from_REF.paf" 2> "$WORK/minimap2.log"
  sort -k1,1 -k3,3n -S 4G "$WORK/from_REF.paf" > "$WORK/from_REF.sorted.paf"
fi
echo "   $(wc -l < "$WORK/from_REF.sorted.paf") alignments"

echo "=== 3. annotation liftover -> GFF3"
# paftools liftover names each lifted interval chr_start_end of its source; the attributes
# are looked up by that key per group (see inputs/gencode_to_groups.sh for why groups)
: > "$WORK/lifted.gff3"
for g in pcgene othergene rest; do
  k8 "$PAFTOOLS" liftover -q 0 -l 1000 -d 1.1 "$WORK/from_REF.sorted.paf" "$GROUPS/$g.bed6" > "$WORK/lifted.$g.bed6"
  awk -F'	' -v OFS='	' 'FNR==NR { typ[$1]=$2; src[$1]=$3; attr[$1]=$4; next }
       ($4 in typ) { print $1, src[$4], typ[$4], $2+1, $3, ".", $6, ".", attr[$4] }' "$GROUPS/$g.map.tsv" "$WORK/lifted.$g.bed6" >> "$WORK/lifted.gff3"
  echo "   $g: lifted $(wc -l < "$WORK/lifted.$g.bed6") of $(wc -l < "$GROUPS/$g.bed6")"
done
awk 'BEGIN{OFS="	"}{print $1,NR}' "$WORK/$SAMPLE.consensus.fa.fai" > "$WORK/order.tsv"
GFF="$DATA/${SAMPLE}_AccessionProcessing/genomic.gff"
{ echo "##gff-version 3"
  awk '{printf("##sequence-region %s 1 %s
",$1,$2)}' "$WORK/$SAMPLE.consensus.fa.fai"
  awk 'BEGIN{FS=OFS="	"} FNR==NR {ord[$1]=$2; next} ($1 in ord){print ord[$1],$0}' "$WORK/order.tsv" "$WORK/lifted.gff3"     | sort -T "$WORK" -k1,1n -k5,5n -S 4G | awk 'BEGIN{OFS="	"}{$1=""; sub(/^	/,""); if(!seen[$0]++){print}}'
} > "$GFF"
echo "   protein-coding genes in genomic.gff: $(awk -F'	' '$3=="gene" && $9 ~ /gene_type=protein_coding/' "$GFF" | wc -l)"

echo "=== 4. chromosome list"
bash "$(dirname "$0")/make_chromosomes_json.sh" "$WORK/$SAMPLE.consensus.fa.fai" "$DATA/${SAMPLE}_AccessionProcessing/${SAMPLE}_Chromosomes.json" "$SAMPLE" "$SAMPLE"
echo "=== done: now run  gammatext-scope prepare-annotation --sample $SAMPLE"
