#!/usr/bin/env bash
# Melanoma v3: same as v2 but with one germline record set, in one order, for both
# samples (see make_vcf_v3.py). Stage 0 builds chr21 for both samples and writes
# it to the Stream *_v3 folders so the breakpoint check can run; PROOF_ONLY=1 stops there.
set -euo pipefail
ROOT="/mnt/d/Cancer sequences"; V2="$ROOT/v2"; V3="$ROOT/v2/v3"; STREAM="/mnt/d/Files/Stream"
BC="$ROOT/breast cancer"; MEL="$ROOT/melanoma"
CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
THREADS=${THREADS:-20}
REF_ENSEMBL="$BC/GRCh38.fa"; REF_MEL="$MEL/GRCh38_chr.fixed.fa"
BED6="$BC/GRCh38.annotation.rich.bed6"; KEYMAP="$BC/grch38.key2name.tsv"; PAFTOOLS="$BC/paftools.js"
mkdir -p "$V3/vcf"
exec > >(tee -a "$V3/prepare_v3.log") 2>&1
echo "=== prepare_v3 start $(date) ==="

SAMPLES=(
  "HomoSapiens_Melanoma_Normal_v3|COLO829BL|$REF_MEL|$V3/vcf/{c}.normal.vcf.gz|-s COLO829 -H A|HomoSapiens_Melanoma_Normal"
  "HomoSapiens_Melanoma_Tumor_v3|COLO829|$REF_MEL|$V3/vcf/{c}.tumour.vcf.gz|-s COLO829 -H A|HomoSapiens_Melanoma_Tumor"
)

echo "=== stage V: per-chromosome VCFs with a shared, ordered germline set ==="
for c in $CHRS; do
  [ -s "$V3/vcf/$c.tumour.vcf.gz" ] && continue
  python3 "$V2/make_vcf_v3.py" "$MEL/tmp_norm_chr/$c.vcf.gz" "$MEL/tmp_som_chr/$c.vcf.gz" "$V3/vcf/$c.normal.vcf" "$V3/vcf/$c.tumour.vcf" | sed "s/^/   chr$c: /"
  for k in normal tumour; do bgzip -f "$V3/vcf/$c.$k.vcf"; tabix -f -p vcf "$V3/vcf/$c.$k.vcf.gz"; done
done

consensus_chr() {  # ref vcfpat opts chr outfa
  local vcf="${2//\{c\}/$4}"
  # shellcheck disable=SC2086
  samtools faidx "$1" "$4" | bcftools consensus $3 "$vcf" > "$5" 2> "$5.log"
  grep -q "does not match" "$5.log" && { echo "   WARNING: REF mismatches in $5.log"; grep -m2 "does not match" "$5.log"; } || true
}

build_chr() {  # spec chr
  IFS='|' read -r name id ref vcfpat opts base <<< "$1"; local c="$2"; local out="$V3/$name"
  mkdir -p "$out/chr" "$STREAM/$name"
  [ -s "$out/chr/chr$c.$id.hap1.fa" ] && return 0
  consensus_chr "$ref" "$vcfpat" "$opts" "$c" "$out/chr/chr$c.$id.hap1.fa"
  grep -v '^>' "$out/chr/chr$c.$id.hap1.fa" > "$STREAM/$name/chr$c.$id.hap1.txt"
}

echo "=== stage 0: chr21 for both samples ==="
for spec in "${SAMPLES[@]}"; do build_chr "$spec" 21; IFS='|' read -r name id _ <<< "$spec"; echo "   $name chr21: $(grep -o 'Applied [0-9]* variants' "$V3/$name/chr/chr21.$id.hap1.fa.log")"; done
[ "${PROOF_ONLY:-0}" = "1" ] && { echo "=== proof only, stopping ==="; exit 0; }

echo "=== stage 1: consensus (all chromosomes) ==="
for spec in "${SAMPLES[@]}"; do
  IFS='|' read -r name id ref vcfpat opts base <<< "$spec"; out="$V3/$name"
  mkdir -p "$STREAM/${name}_AccessionProcessing" "$STREAM/${name}_GeneInfo"
  for c in $CHRS; do build_chr "$spec" "$c"; done
  files=(); for c in $CHRS; do files+=("$out/chr/chr$c.$id.hap1.fa"); done
  cat "${files[@]}" > "$out/$name.consensus.fa"; samtools faidx "$out/$name.consensus.fa"
  echo "   $name: $(ls "$STREAM/$name" | wc -l) chromosome files; $(awk '{s+=$2}END{print s}' "$out/$name.consensus.fa.fai") bp ($(date +%H:%M))"
done

echo "=== stage 2: minimap2 asm5 (reference as query) ==="
align() { local name="$1"; local out="$V3/$name"
  [ -s "$out/${name}_from_REF.cg.sorted.paf" ] && return 0
  minimap2 -t "$THREADS" -x asm5 -c --secondary=no "$out/$name.consensus.fa" "$REF_ENSEMBL" > "$out/${name}_from_REF.cg.paf" 2> "$out/minimap2.log"
  sort -k1,1 -k3,3n -S 4G "$out/${name}_from_REF.cg.paf" > "$out/${name}_from_REF.cg.sorted.paf"
  echo "   $name: $(wc -l < "$out/${name}_from_REF.cg.sorted.paf") alignments ($(date +%H:%M))"; }
align HomoSapiens_Melanoma_Normal_v3 & p1=$!; align HomoSapiens_Melanoma_Tumor_v3 & p2=$!; wait $p1 $p2

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
  IFS='|' read -r name id ref vcfpat opts base <<< "$spec"; out="$V3/$name"; acc="$STREAM/${name}_AccessionProcessing"
  [ -s "$acc/genomic.gff" ] && continue
  k8 "$PAFTOOLS" liftover -q 0 -l 1000 -d 1.1 "$out/${name}_from_REF.cg.sorted.paf" "$BED6" > "$out/$name.annotation.rich.bed6"
  lifted_bed_to_gff "$out/$name.annotation.rich.bed6" "$KEYMAP" "$out/$name.annotation.gff3"
  awk 'BEGIN{OFS="\t"}{print $1,NR}' "$out/$name.consensus.fa.fai" > "$out/order.tsv"
  { echo "##gff-version 3"; awk '{printf("##sequence-region %s 1 %s\n",$1,$2)}' "$out/$name.consensus.fa.fai"
    awk 'BEGIN{FS=OFS="\t"} FNR==NR {ord[$1]=$2; next} /^#/ {next} ($1 in ord){print ord[$1],$0}' "$out/order.tsv" "$out/$name.annotation.gff3" \
      | sort -T "$out" -k1,1n -k5,5n -S 4G | awk 'BEGIN{OFS="\t"}{$1=""; sub(/^\t/,""); if(!seen[$0]++){print}}'
  } > "$acc/genomic.gff"
  cp "$STREAM/${base}_AccessionProcessing/${base}_Chromosomes.json" "$acc/${name}_Chromosomes.json"
  echo "   $name: $(wc -l < "$out/$name.annotation.rich.bed6") lifted features; protein-coding genes: $(awk -F'\t' '$3=="gene" && $9 ~ /gene_type=protein_coding/' "$acc/genomic.gff" | wc -l)"
done
echo "=== prepare_v3 done $(date) ==="
