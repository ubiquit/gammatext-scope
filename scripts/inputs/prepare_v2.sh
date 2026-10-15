#!/usr/bin/env bash
# gammatext-scope article-9 inputs, version 2 (2026-10-15).
# Builds, for all four human samples, a pseudo-haploid consensus genome from the
# existing VCFs with ONE explicit rule (bcftools consensus -H A: the ALT allele of
# every genotype), and lifts the GENCODE v44 annotation onto each consensus with
# the breast-cancer recipe (minimap2 asm5, paftools.js liftover, correct direction:
# the reference is the query). Writes the sequences and annotation straight into
# the GenomeProcessor/gammatext-scope data layout under new *_v2 sample names.
#
# Stage 0 is a proof: chromosome 21 of every sample is rebuilt and compared with
# the existing (v1) string. BC normal, BC tumour and melanoma tumour must be
# identical (same rule as before); melanoma normal must differ (that is the fix).
# Any other outcome aborts the run.
set -euo pipefail
ROOT="/mnt/d/Cancer sequences"; V2="$ROOT/v2"; STREAM="/mnt/d/Files/Stream"
BC="$ROOT/breast cancer"; MEL="$ROOT/melanoma"
CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
THREADS=${THREADS:-20}
REF_ENSEMBL="$BC/GRCh38.fa"                 # Ensembl GRCh38 primary assembly, names 1..22 X Y MT
REF_MEL="$MEL/GRCh38_chr.fixed.fa"          # reference the melanoma VCFs were normalised against (same primary sequence)
BED6="$BC/GRCh38.annotation.rich.bed6"      # GENCODE v44 (Ensembl 110) as BED6 with TYPE/SOURCE/attributes in the name column
KEYMAP="$BC/grch38.key2name.tsv"
PAFTOOLS="$BC/paftools.js"                  # 2.30-r1287
mkdir -p "$V2"
exec > >(tee -a "$V2/prepare_v2.log") 2>&1
echo "=== prepare_v2 start $(date) === threads=$THREADS"
for t in bcftools samtools minimap2 k8; do printf "%s: %s\n" "$t" "$($t --version 2>&1 | head -1)"; done

# name | sequence id | reference | per-chromosome VCF pattern ({c}) | bcftools consensus options | v1 sample name
SAMPLES=(
  "HomoSapiens_BreastCancer_Normal_v2|HCC1395BL|$REF_ENSEMBL|$BC/HCC1395BL_germline.norm.vcf.gz|-s sample -H A|HomoSapiens_BreastCancer_Normal"
  "HomoSapiens_BreastCancer_Tumor_v2|HCC1395|$REF_ENSEMBL|$BC/tumor_union.vcf.gz||HomoSapiens_BreastCancer_Tumor"
  "HomoSapiens_Melanoma_Normal_v2|COLO829BL|$REF_MEL|$MEL/tmp_norm_chr/{c}.vcf.gz|-s COLO829 -H A|HomoSapiens_Melanoma_Normal"
  "HomoSapiens_Melanoma_Tumor_v2|COLO829|$REF_MEL|$MEL/tmp_union_chr/{c}.union.vcf.gz|-s COLO829 -H A|HomoSapiens_Melanoma_Tumor"
)

consensus_chr() {  # ref vcfpat opts chr outfa
  local ref="$1" vcfpat="$2" opts="$3" c="$4" out="$5"
  local vcf="${vcfpat//\{c\}/$c}"
  # shellcheck disable=SC2086
  samtools faidx "$ref" "$c" | bcftools consensus $opts "$vcf" > "$out" 2> "$out.log"
  grep -c "does not match" "$out.log" >/dev/null 2>&1 && { echo "   WARNING: REF mismatches in $out.log"; grep -m2 "does not match" "$out.log"; } || true
}

# ---------------------------------------------------------------- stage 0: proof on chr21
echo "=== stage 0: chr21 proof ==="
mkdir -p "$V2/proof"
for spec in "${SAMPLES[@]}"; do
  IFS='|' read -r name id ref vcfpat opts base <<< "$spec"
  consensus_chr "$ref" "$vcfpat" "$opts" 21 "$V2/proof/chr21.$id.fa"
  new=$(grep -v '^>' "$V2/proof/chr21.$id.fa" | md5sum | cut -c1-32)
  old=$(md5sum "$STREAM/$base/chr21.$id.hap1.txt" | cut -c1-32)
  if [ "$new" = "$old" ]; then verdict="IDENTICAL to v1"; else verdict="DIFFERENT from v1"; fi
  echo "   $name chr21: $verdict  (applied: $(grep -o 'Applied [0-9]* variants' "$V2/proof/chr21.$id.fa.log" || echo '?'))"
  case "$name" in
    *Melanoma_Normal*) [ "$new" != "$old" ] || { echo "ABORT: melanoma normal chr21 did not change — the fix had no effect"; exit 2; } ;;
    *)                 [ "$new" = "$old" ]  || { echo "ABORT: $name chr21 differs from v1 although the rule should reproduce it"; exit 2; } ;;
  esac
done
[ "${PROOF_ONLY:-0}" = "1" ] && { echo "=== proof only, stopping ==="; exit 0; }

# ---------------------------------------------------------------- stage 1: consensus genomes
echo "=== stage 1: consensus (all chromosomes) ==="
for spec in "${SAMPLES[@]}"; do
  IFS='|' read -r name id ref vcfpat opts base <<< "$spec"
  out="$V2/$name"; mkdir -p "$out/chr" "$STREAM/$name" "$STREAM/${name}_AccessionProcessing" "$STREAM/${name}_GeneInfo"
  echo "--- $name ($(date +%H:%M))"
  for c in $CHRS; do
    [ -s "$out/chr/chr$c.$id.hap1.fa" ] && continue
    consensus_chr "$ref" "$vcfpat" "$opts" "$c" "$out/chr/chr$c.$id.hap1.fa"
    grep -v '^>' "$out/chr/chr$c.$id.hap1.fa" > "$STREAM/$name/chr$c.$id.hap1.txt"
  done
  files=(); for c in $CHRS; do files+=("$out/chr/chr$c.$id.hap1.fa"); done
  cat "${files[@]}" > "$out/$name.consensus.fa"
  samtools faidx "$out/$name.consensus.fa"
  echo "   $(ls "$STREAM/$name" | wc -l) chromosome files; total $(awk '{s+=$2}END{print s}' "$out/$name.consensus.fa.fai") bp"
done

# ---------------------------------------------------------------- stage 2: alignment consensus <- reference (two at a time)
echo "=== stage 2: minimap2 asm5 (reference as query) ==="
align() { local name="$1"; local out="$V2/$name"
  [ -s "$out/${name}_from_REF.cg.sorted.paf" ] && return 0
  minimap2 -t "$THREADS" -x asm5 -c --secondary=no "$out/$name.consensus.fa" "$REF_ENSEMBL" > "$out/${name}_from_REF.cg.paf" 2> "$out/minimap2.log"
  sort -k1,1 -k3,3n -S 4G "$out/${name}_from_REF.cg.paf" > "$out/${name}_from_REF.cg.sorted.paf"
  echo "   $name: $(wc -l < "$out/${name}_from_REF.cg.sorted.paf") alignments ($(date +%H:%M))"; }
align HomoSapiens_BreastCancer_Normal_v2 & p1=$!; align HomoSapiens_BreastCancer_Tumor_v2 & p2=$!; wait $p1 $p2   # bare 'wait' would also wait for the tee process substitution (bash 5.1)
align HomoSapiens_Melanoma_Normal_v2     & p1=$!; align HomoSapiens_Melanoma_Tumor_v2     & p2=$!; wait $p1 $p2

# ---------------------------------------------------------------- stage 3: annotation liftover (breast-cancer recipe)
echo "=== stage 3: paftools liftover + GFF3 ==="
lifted_bed_to_gff () {   inbed="$1"; keymap="$2"; outgff="$3";   printf "##gff-version 3\n" > "$outgff";   awk -v FS="\t" -v OFS="\t" '
    FNR==NR { map[$1]=$2; next }
    {
      tchr=$1; tstart=$2; tend=$3; key=$4; strand=$6;
      name=map[key]; if(name=="") name=".";
      typ="feature"; src="liftover"; attrs=name;
      if (match(name, /(^|;)TYPE=([^;]+)/, m)) { typ=m[2] }
      if (match(name, /(^|;)SOURCE=([^;]+)/, m)) { src=m[2] }
      gsub(/(^|;)TYPE=[^;]+;?/, "", attrs);
      gsub(/(^|;)SOURCE=[^;]+;?/, "", attrs);
      sub(/^;/, "", attrs); sub(/;$/, "", attrs); if (attrs=="") attrs=".";
      print tchr, src, typ, tstart+1, tend, ".", strand, ".", attrs
    }' "$keymap" "$inbed" >> "$outgff"; }

for spec in "${SAMPLES[@]}"; do
  IFS='|' read -r name id ref vcfpat opts base <<< "$spec"
  out="$V2/$name"; acc="$STREAM/${name}_AccessionProcessing"
  [ -s "$acc/genomic.gff" ] && { echo "   $name: genomic.gff exists, skipping"; continue; }
  k8 "$PAFTOOLS" liftover -q 0 -l 1000 -d 1.1 "$out/${name}_from_REF.cg.sorted.paf" "$BED6" > "$out/$name.annotation.rich.bed6"
  lifted_bed_to_gff "$out/$name.annotation.rich.bed6" "$KEYMAP" "$out/$name.annotation.gff3"
  awk 'BEGIN{OFS="\t"}{print $1,NR}' "$out/$name.consensus.fa.fai" > "$out/order.tsv"
  { echo "##gff-version 3"
    awk '{printf("##sequence-region %s 1 %s\n",$1,$2)}' "$out/$name.consensus.fa.fai"
    awk 'BEGIN{FS=OFS="\t"} FNR==NR {ord[$1]=$2; next} /^#/ {next} ($1 in ord){print ord[$1],$0}' "$out/order.tsv" "$out/$name.annotation.gff3" \
      | sort -T "$out" -k1,1n -k5,5n -S 4G | awk 'BEGIN{OFS="\t"}{$1=""; sub(/^\t/,""); if(!seen[$0]++){print}}'
  } > "$acc/genomic.gff"
  cp "$STREAM/${base}_AccessionProcessing/${base}_Chromosomes.json" "$acc/${name}_Chromosomes.json"
  echo "   $name: $(wc -l < "$out/$name.annotation.rich.bed6") lifted features of $(wc -l < "$BED6"); protein-coding genes in genomic.gff: $(awk -F'\t' '$3=="gene" && $9 ~ /gene_type=protein_coding/' "$acc/genomic.gff" | wc -l)"
done
echo "=== prepare_v2 done $(date) ==="
