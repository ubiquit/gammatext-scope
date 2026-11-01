# Building the input genomes

How the four sample folders under the data root (see README, *Data layout*)
are built from public data. The scripts in `inputs/` do the generic steps;
the variant preparation is a handful of bcftools commands listed here.

Tools used: bcftools 1.13, samtools 1.13, minimap2 2.24 with its
`paftools.js` (run with `k8`), Python 3.

## Sources

| | breast cancer (HCC1395 / HCC1395BL) | melanoma (COLO829 / COLO829BL) |
|---|---|---|
| reference | GRCh38 primary assembly from Ensembl, sequences `1..22 X Y MT` (`GRCh38.fa`) | same sequence; the ONT VCFs name the contigs `chr1..chrY, chrM`, see the rename step |
| normal | NCBI assembly **GCA_021234545.1** (HCC1395BL_v1.0) | ONT open data, `COLO829_normal_germline.vcf.gz` (Clair3 1.0.4, sup basecalling) |
| tumour | SEQC2 `high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` and `high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz` (Fang et al. 2021) | ONT open data, `COLO829_somatic_snv.vcf.gz` and `COLO829_somatic_indels.vcf.gz` (ClairS 0.1.7) |
| annotation | GENCODE release 44, `gencode.v44.annotation.gff3.gz` | same |

`rename.chr.map` is a two-column file `chr1 1` … `chrX X`, `chrY Y`, `chrM MT`.

## 1. Variants

**Breast-cancer normal** — germline variants are called from the assembly
aligned to the reference:

```bash
minimap2 -cx asm5 --cs -t 20 GRCh38.fa HCC1395BL_assembly.fa > BLvsGRCh38.paf
sort -k6,6 -k8,8n BLvsGRCh38.paf | k8 paftools.js call -f GRCh38.fa -s sample - > HCC1395BL_germline.raw.vcf
bcftools sort -Oz -o HCC1395BL_germline.raw.sorted.vcf.gz HCC1395BL_germline.raw.vcf && tabix HCC1395BL_germline.raw.sorted.vcf.gz
bcftools norm -f GRCh38.fa -Oz -o HCC1395BL_germline.norm.vcf.gz HCC1395BL_germline.raw.sorted.vcf.gz && tabix HCC1395BL_germline.norm.vcf.gz
```

**Breast-cancer tumour** — the germline records plus the SEQC2 somatic
records, as a sites-only union:

```bash
for k in sSNV sINDEL; do
  bcftools annotate --rename-chrs rename.chr.map -Oz -o highconf_${k}_renamed.vcf.gz high-confidence_${k}_in_HC_regions_v1.2.1.vcf.gz && tabix highconf_${k}_renamed.vcf.gz
done
bcftools concat -a highconf_sSNV_renamed.vcf.gz highconf_sINDEL_renamed.vcf.gz | bcftools norm -f GRCh38.fa -Oz -o tumor_somatic.norm.vcf.gz && tabix tumor_somatic.norm.vcf.gz
bcftools merge -m none -Oz -o tumor_union.vcf.gz HCC1395BL_germline.norm.vcf.gz tumor_somatic.norm.vcf.gz && tabix tumor_union.vcf.gz
```

**Melanoma, both samples** — normalise, split per chromosome, rename the
contigs, then give both samples one germline record set in one order with the
somatic records merged into the tumour copy (`inputs/make_vcf_v3.py`; this is
what makes `bcftools consensus` take the same apply/skip decisions in both
samples):

```bash
bcftools norm -f GRCh38.fa -m -any -Oz -o normal.germline.norm.vcf.gz COLO829_normal_germline.vcf.gz && tabix normal.germline.norm.vcf.gz
bcftools concat -a COLO829_somatic_snv.vcf.gz COLO829_somatic_indels.vcf.gz | bcftools view -f PASS \
  | bcftools norm -f GRCh38.fa -m -any | bcftools sort -Oz -o tumor.somatic.PASS.norm.sorted.vcf.gz && tabix tumor.somatic.PASS.norm.sorted.vcf.gz
mkdir -p norm_chr som_chr vcf
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
  bcftools view -r chr$c normal.germline.norm.vcf.gz          | bcftools annotate --rename-chrs rename.chr.map -Oz -o norm_chr/$c.vcf.gz && tabix norm_chr/$c.vcf.gz
  bcftools view -r chr$c tumor.somatic.PASS.norm.sorted.vcf.gz | bcftools annotate --rename-chrs rename.chr.map -Oz -o som_chr/$c.vcf.gz  && tabix som_chr/$c.vcf.gz
  python3 inputs/make_vcf_v3.py norm_chr/$c.vcf.gz som_chr/$c.vcf.gz vcf/$c.normal.vcf vcf/$c.tumour.vcf
  for k in normal tumour; do bgzip -f vcf/$c.$k.vcf && tabix -f -p vcf vcf/$c.$k.vcf.gz; done
done
```

## 2. Annotation as BED6

```bash
inputs/gencode_to_groups.sh gencode.v44.annotation.gff3.gz groups
```

writes three BED6 files with a key map each under `groups/`: protein-coding
genes, other genes, everything else. `paftools liftover` names every lifted
interval by its source coordinates, so the attributes are restored by
coordinates afterwards; lifting the genes as their own groups keeps a gene
from being overwritten by a transcript of the same span.

## 3. The four samples

`inputs/build_sample.sh` builds one sample: the consensus of every chromosome
(`bcftools consensus` with the ALT allele of every genotype), the alignment of
the reference to the consensus, the liftover of the annotation and the
chromosome list. Set `GAMMATEXT_DATA` to the data root first.

```bash
inputs/build_sample.sh HomoSapiens_BreastCancer_Normal_v2 HCC1395BL GRCh38.fa HCC1395BL_germline.norm.vcf.gz "-s sample -H A" groups paftools.js
inputs/build_sample.sh HomoSapiens_BreastCancer_Tumor_v2  HCC1395   GRCh38.fa tumor_union.vcf.gz              ""              groups paftools.js
inputs/build_sample.sh HomoSapiens_Melanoma_Normal_v3     COLO829BL GRCh38.fa "vcf/{c}.normal.vcf.gz"         "-s COLO829 -H A" groups paftools.js
inputs/build_sample.sh HomoSapiens_Melanoma_Tumor_v3      COLO829   GRCh38.fa "vcf/{c}.tumour.vcf.gz"         "-s COLO829 -H A" groups paftools.js
```

The breast-cancer tumour union carries no genotypes, so `bcftools consensus`
applies the ALT allele of every record without `-H`.

Then make the gene mask of each pair identical. The annotation is lifted onto
each consensus separately, and paftools can reject a gene in one sample only
(an indel inside the gene breaks its length check); such a gene would be
masked in one filtered sequence and not in the other, and every window after
it on that chromosome would compare different DNA. The script removes the
protein-coding `gene` lines present in only one file of the pair:

```bash
python inputs/pair_masks.py data/HomoSapiens_BreastCancer_Normal_v2_AccessionProcessing/genomic.gff data/HomoSapiens_BreastCancer_Tumor_v2_AccessionProcessing/genomic.gff
python inputs/pair_masks.py data/HomoSapiens_Melanoma_Normal_v3_AccessionProcessing/genomic.gff     data/HomoSapiens_Melanoma_Tumor_v3_AccessionProcessing/genomic.gff
```

For dataset v2 this removes nothing from the breast-cancer pair (19,743
protein-coding genes in both) and one gene, ABR on chr17, from the melanoma
normal (19,290 in both afterwards).

## 4. The pipeline

```bash
for s in HomoSapiens_BreastCancer_Normal_v2 HomoSapiens_BreastCancer_Tumor_v2 HomoSapiens_Melanoma_Normal_v3 HomoSapiens_Melanoma_Tumor_v3; do
  gammatext-scope prepare-annotation --sample $s
  for m in Intergenic Gene WholeGenome; do gammatext-scope classify --sample $s --mode $m; done
done
for m in Intergenic Gene WholeGenome; do
  gammatext-scope density       --normal HomoSapiens_BreastCancer_Normal_v2 --tumor HomoSapiens_BreastCancer_Tumor_v2 --pair BreastCancer_v2_Normal_vs_Tumor --mode $m
  gammatext-scope alpha-compare --normal HomoSapiens_BreastCancer_Normal_v2 --tumor HomoSapiens_BreastCancer_Tumor_v2 --pair BreastCancer_v2 --mode $m
  gammatext-scope decay-dfa     --pair BreastCancer_v2 --mode $m
  gammatext-scope density       --normal HomoSapiens_Melanoma_Normal_v3 --tumor HomoSapiens_Melanoma_Tumor_v3 --pair Melanoma_v3_Normal_vs_Tumor --mode $m
  gammatext-scope alpha-compare --normal HomoSapiens_Melanoma_Normal_v3 --tumor HomoSapiens_Melanoma_Tumor_v3 --pair Melanoma_v3 --mode $m
  gammatext-scope decay-dfa     --pair Melanoma_v3 --mode $m
done
gammatext-scope alpha-compare --normal HomoSapiens_BreastCancer_Normal_v2 --tumor HomoSapiens_BreastCancer_Tumor_v2 --pair BreastCancer_v2 --mode WholeGenome --merged
gammatext-scope decay-dfa     --pair BreastCancer_v2 --mode WholeGenome --merged
gammatext-scope alpha-compare --normal HomoSapiens_Melanoma_Normal_v3 --tumor HomoSapiens_Melanoma_Tumor_v3 --pair Melanoma_v3 --mode WholeGenome --merged
gammatext-scope decay-dfa     --pair Melanoma_v3 --mode WholeGenome --merged
```

## 5. Checks

* `density` and `alpha-compare` print the number of alignment breakpoints per
  chromosome; on chr21 it equals the somatic indel count of the pair (breast
  cancer 21, melanoma 23).
* At heterozygous and homozygous germline sites both consensus genomes carry
  the ALT allele; at somatic sites only the tumour does.
* `build_sample.sh` prints the number of protein-coding genes in each lifted
  annotation; GENCODE v44 has 20,031 on the primary chromosomes, and the
  liftover (`-l 1000 -d 1.1`) loses 288 of them in the breast-cancer samples
  and 740 in the melanoma samples. `pair_masks.py` reports what it removes.
