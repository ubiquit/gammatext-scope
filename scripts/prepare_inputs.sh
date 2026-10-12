#!/usr/bin/env bash
# Reconstructs the per-chromosome input sequences and lifted annotations used by
# the article from public data. Reconstructed from the bcftools command lines
# stamped in the intermediate VCF headers and the artefacts left on disk;
# lines marked TODO(verify) are steps whose exact flags were not recoverable
# and must be confirmed against the original run before this script is relied on.
#
# Tools (versions used originally): bcftools 1.13, minimap2/paftools.js 2.30-r1287,
# k8 (b52d6a8), CrossMap, samtools. Run under Linux/WSL.
set -euo pipefail

REF=GRCh38.fa                      # Ensembl GRCh38 primary assembly, headers ">1 dna:chromosome …"
GENCODE=gencode.v44.annotation.gff3 # GENCODE v44 (Ensembl 110), chr-prefixed
MAP=rename.chr.map                 # "chr1 1" … "chrM MT"

# ---------------------------------------------------------------- breast cancer
# Tumour: SEQC2 high-confidence somatic calls v1.2.1 (Fang et al. 2021)
bcftools annotate --rename-chrs "$MAP" -Oz -o highconf_sSNV_renamed.vcf.gz   high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
bcftools annotate --rename-chrs "$MAP" -Oz -o highconf_sINDEL_renamed.vcf.gz high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
bcftools concat -a highconf_sSNV_renamed.vcf.gz highconf_sINDEL_renamed.vcf.gz -Oz -o tumor_somatic.merged.vcf.gz
bcftools norm -f "$REF" -Oz -o HCC1395_tumor_somatic.norm.vcf.gz tumor_somatic.merged.vcf.gz
bcftools index -t HCC1395_tumor_somatic.norm.vcf.gz

# Normal: differences between the NCBI de novo assembly of HCC1395BL and GRCh38
# (NCBI datasets GCA_021234545.1, HCC1395BL_v1.0)
minimap2 -cx asm5 --cs -t 8 "$REF" GCA_021234545.1_HCC1395BL_v1.0_genomic.fna > BLvsGRCh38.paf   # TODO(verify) preset
sort -k6,6 -k8,8n BLvsGRCh38.paf | k8 paftools.js call -f "$REF" - > HCC1395BL_germline.raw.vcf  # TODO(verify) flags
bcftools sort -Oz -o HCC1395BL_germline.raw.sorted.vcf.gz HCC1395BL_germline.raw.vcf
bcftools norm -f "$REF" -Oz -o HCC1395BL_germline.norm.vcf.gz HCC1395BL_germline.raw.sorted.vcf.gz
bcftools index -t HCC1395BL_germline.norm.vcf.gz

# Tumour union: germline calls minus somatic positions, plus somatic sites
bcftools query -f '%CHROM\t%POS\n' HCC1395_tumor_somatic.norm.vcf.gz > somatic.pos
bcftools view -T ^somatic.pos -Oz -o germline.noSom.vcf.gz HCC1395BL_germline.norm.vcf.gz
bcftools view -G -Oz -o germline.noSom.sites.vcf.gz germline.noSom.vcf.gz
bcftools view -G -Oz -o tumor_somatic.sites.vcf.gz HCC1395_tumor_somatic.norm.vcf.gz
bcftools concat -a -Oz -o tumor_union.vcf.gz tumor_somatic.sites.vcf.gz germline.noSom.sites.vcf.gz
bcftools sort -Oz -o tumor_union.sorted.vcf.gz tumor_union.vcf.gz && bcftools index -t tumor_union.sorted.vcf.gz

# Consensus sequences, one file per chromosome, FASTA header stripped
for c in $(seq 1 22) X Y MT; do
  bcftools view -r "$c" -Oz -o tmp.normal.$c.vcf.gz HCC1395BL_germline.norm.vcf.gz && bcftools index -t tmp.normal.$c.vcf.gz
  samtools faidx "$REF" "$c" | bcftools consensus -H 1 tmp.normal.$c.vcf.gz > chr$c.HCC1395BL.hap1.fa   # TODO(verify) -H
  bcftools view -r "$c" -Oz -o tmp.tumor.$c.vcf.gz tumor_union.sorted.vcf.gz && bcftools index -t tmp.tumor.$c.vcf.gz
  samtools faidx "$REF" "$c" | bcftools consensus tmp.tumor.$c.vcf.gz > chr$c.HCC1395.hap1.fa              # TODO(verify)
  tail -n +2 chr$c.HCC1395BL.hap1.fa > HomoSapiens_BreastCancer_Normal/chr$c.HCC1395BL.hap1.txt
  tail -n +2 chr$c.HCC1395.hap1.fa   > HomoSapiens_BreastCancer_Tumor/chr$c.HCC1395.hap1.txt
done

# Annotation liftover onto each consensus (GENCODE v44 → BED6 → paftools liftover → GFF3)
# TODO(verify): BED6 conversion, minimap2 alignment of consensus vs reference (REF_to_NORMAL.cg.paf), paftools liftover, BED→GFF3

# ---------------------------------------------------------------- melanoma (ONT open data COLO829 / COLO829BL)
bcftools norm -f "$REF" -m -any -Oz -o normal.germline.norm.vcf.gz COLO829_normal_germline.vcf.gz          # Clair3 1.0.4
bcftools concat -a COLO829_somatic_snv.vcf.gz COLO829_somatic_indels.vcf.gz | bcftools view -f PASS | \
  bcftools norm -f "$REF" -m -any | bcftools sort -Oz -o tumor.somatic.PASS.norm.sorted.vcf.gz             # ClairS 0.1.7
bcftools query -f '%CHROM\t%POS0\t%END\n' tumor.somatic.PASS.norm.sorted.vcf.gz > somatic.spans.bed
bcftools view -T ^somatic.spans.bed -Oz -o normal.germline.noSomOverlapSpan.vcf.gz normal.germline.norm.vcf.gz
for c in $(seq 1 22) X Y MT; do
  bcftools view -r chr$c normal.germline.norm.vcf.gz | bcftools annotate --rename-chrs rename.toFASTA.map -Oz -o tmp_norm_chr/$c.vcf.gz
  bcftools view -r chr$c tumor.somatic.PASS.norm.sorted.vcf.gz | bcftools annotate --rename-chrs rename.toFASTA.map -Oz -o tmp_som_chr/$c.vcf.gz
  bcftools concat -a -Oz -o tmp_union_chr/$c.vcf.gz tmp_union_chr/$c.germ.vcf.gz tmp_union_chr/$c.som.vcf.gz
  samtools faidx "$REF" "$c" | bcftools consensus -c chains_normal/$c.chain -H 1 tmp_norm_chr/$c.vcf.gz > chr$c.COLO829BL.hap1.fa   # TODO(verify) -H
  samtools faidx "$REF" "$c" | bcftools consensus -c chains_tumor/$c.chain tmp_union_chr/$c.vcf.gz       > chr$c.COLO829.hap1.fa     # TODO(verify)
done
# Annotation liftover: CrossMap.py gff <chain> gencode.v44.noprefix.gff3 <out.gff3> per sample   # TODO(verify)
