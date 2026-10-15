# Input preparation — the scripts that were actually run

These are the scripts, unchanged, that produced the sequences and annotations
under `expected/v2` on 2026-10-15 (Linux side under WSL Ubuntu 22.04: bcftools
1.13, samtools 1.13, minimap2 2.24, paftools.js 2.30-r1287 on k8; Windows side:
gammatext-scope). Paths are the authors' and need adapting.

| script | what it does |
|---|---|
| `run_v2.sh` | driver: backup of v1 → `prepare_v2.sh` → `prepare-annotation`, `classify`, experiments → `verify_v2.py` |
| `prepare_v2.sh` | consensus genomes for all four samples with one rule (`bcftools consensus -H A`), annotation lifted per sample; chr21 proof against v1 first |
| `make_vcf_v3.py`, `prepare_v3.sh`, `run_v3.sh` | melanoma only: one ordered germline record set for both samples (records overlapping a somatic span removed), somatic merged into the tumour copy; this is the melanoma dataset used (`*_v3`) |
| `verify_v2.py` | allele-application check at heterozygous, homozygous and somatic sites; breast-cancer identity with v1 |

## Construction rule (dataset v2)

1. **Variants.** Breast cancer: normal = `paftools.js call` of the NCBI assembly
   GCA_021234545.1 (HCC1395BL_v1.0) against GRCh38; tumour = normal minus the
   SEQC2 somatic positions plus the SEQC2 high-confidence sSNV + sINDEL v1.2.1
   (sites-only union). Melanoma: ONT open data, Clair3 1.0.4 germline and
   ClairS 0.1.7 somatic (PASS), `bcftools norm -m -any`; normal and tumour
   share one germline record set in one order (`make_vcf_v3.py`).
2. **Consensus.** `samtools faidx GRCh38 <chr> | bcftools consensus -H A <vcf>`
   per chromosome — the ALT allele of every genotype, on both samples
   (breast-cancer tumour union has no genotypes; bcftools applies ALT).
   FASTA header stripped → `chr<N>.<id>.hap1.txt`.
3. **Annotation.** GENCODE v44 (Ensembl 110) as BED6 →
   `minimap2 -x asm5 -c --secondary=no <consensus> GRCh38` (reference as the
   query) → `paftools.js liftover -q 0 -l 1000 -d 1.1` → GFF3, per sample.
4. **Checks.** Normal–tumour alignment breakpoints on chr21 equal the somatic
   indel count (breast cancer 21, melanoma 23); heterozygous and homozygous
   germline alleles present on both sides, somatic alleles only in the tumour.

## What was wrong with v1 (the submitted manuscript's inputs)

* Melanoma normal built with `-H 1` (reference allele at unphased
  heterozygous sites) while the tumour carried the ALT alleles — ~97 % of the
  tumour–normal difference was germline heterozygosity.
* Melanoma tumour union ordered multi-allelic records opposite to the normal,
  so `bcftools consensus` applied different alleles at ~1,200 indel sites per
  chromosome 21.
* Both breast-cancer samples used one annotation lifted (with the alignment
  direction inverted) onto an earlier tumour FASTA without the germline
  insertions, i.e. in a coordinate system up to ~100 kb off the analysed
  strings; gene boundaries were therefore cut at non-homologous positions in
  normal vs tumour.
