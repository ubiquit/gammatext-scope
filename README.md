# gammatext-scope

Stop-codon γtext comparison of a tumour genome with its matched normal genome —
the computational pipeline behind the **DNA-scope** described in

> Glattre H.R., Glattre E., Moe L. *Understanding Carcinogenesis Differently — Our
> Stop-Codon World / A General Cancer Causation Formula* (manuscript).

Non-protein-coding DNA is read as γtext: strings of nucleotides delimited by the
three stop codons TAA, TAG and TGA are γwords. The pipeline aligns a tumour
chromosome to its normal counterpart, counts stop codons in windows of 7,000
γwords, classifies each window by the Zipf/DFA behaviour of its words
(*zipf*, *alfa* or unqualified), and submits the series of qualification-decay
markers to detrended fluctuation analysis (DFA).

This repository is a self-contained lift of the article's code path from the
authors' internal GenomeProcessor platform (commit `8b3c1a8`). The analysis
code is kept verbatim where it affects results; only I/O locations, chart
labels and dead output were changed (see *Differences from the legacy
implementation*).

## Pipeline

| command | input | output | used for |
|---|---|---|---|
| `prepare-annotation` | `<S>_AccessionProcessing/genomic.gff` + `<S>_Chromosomes.json` | `<S>_GeneInfo/<S>_Annotation_Converted.txt` | protein-coding gene mask |
| `classify` | per-chromosome sequence + gene mask | `Experiment_<S>_<mode>_7000/WindowWidth1/chr*WW1_window_positions_indexLength.txt` | window qualification (*zipf* / *alfa* / unqualified) |
| `density` | normal + tumour sequences | `StopCodonDensity_<mode>/*_perChromosome.tsv` … | Fig. 1 (per-chromosome stop-codon difference) |
| `alpha-compare` | sequences + both samples' classifications | `StopCodonAlphaComparison_<mode>/StopCodonAlpha_<pair>_chr*.tsv`, `_AllChromosomes.tsv` | DNA-scopes (Fig. 2–3, Appendices 1–2), Table 3 window counts |
| `decay-dfa` | `_AllChromosomes.tsv` | `DecayGapDFA_<mode>/DecayGapDFA_<pair>_Summary.txt` | Tables 2–3 (DFA α of decay-marker gaps) |

`<mode>` is `Intergenic` (npc-DNA), `Gene` (pc-DNA) or `WholeGenome` (the
unmasked chromosome). For `alpha-compare` and `decay-dfa`, `--mode WholeGenome`
requires whole-genome window classifications for both samples; `--merged`
selects the legacy behaviour instead — the union of the Intergenic and Gene
windows in genomic order (outputs suffixed `WholeGenomeMerged`), which depends
on the gene annotation.

## Build and run

Requires the .NET 8 SDK.

```bash
dotnet build -c Release
dotnet test
src/GammaTextScope/bin/Release/net8.0/gammatext-scope --help
```

Input and output roots default to the authors' layout (`D:\Files\Stream`,
`D:\Files\Harvests`) and are overridden with `--data-root` / `--out-root`
(or `GAMMATEXT_DATA` / `GAMMATEXT_OUT`). `--chromosomes 21,22` restricts any
command to a subset of chromosomes.

### Data layout

```
<data-root>/
  HomoSapiens_BreastCancer_Normal/            chr1.HCC1395BL.hap1.txt … chrY.…   (sequence only, no FASTA header)
  HomoSapiens_BreastCancer_Normal_GeneInfo/   HomoSapiens_BreastCancer_Normal_Annotation_Converted.txt
  HomoSapiens_BreastCancer_Normal_AccessionProcessing/
                                              genomic.gff, HomoSapiens_BreastCancer_Normal_Chromosomes.json
  HomoSapiens_BreastCancer_Tumor/ …           (same for the tumour and for the melanoma pair)
```

Only `gene` features with `gene_type=protein_coding` are read from the GFF.

## Reproducing the article

`expected/` holds two reference datasets — tables, DFA summaries and the
per-window classifications of all four samples — so the experiments can be
re-run and diffed without re-classifying:

* `expected/v1` — the outputs the submitted manuscript (October 2026) was
  built from. Kept for the record; its inputs have known construction defects
  (see *Input provenance*).
* `expected/v2` — the regenerated dataset: breast-cancer samples `*_v2`,
  melanoma samples `*_v3`, pairs `BreastCancer_v2*` / `Melanoma_v3*`. This is
  the dataset the revised manuscript uses.

```bash
# Windows PowerShell / bash: see scripts/reproduce_article9.sh
gammatext-scope density       --normal HomoSapiens_BreastCancer_Normal --tumor HomoSapiens_BreastCancer_Tumor --pair BreastCancer_Normal_vs_Tumor --mode Intergenic
gammatext-scope alpha-compare --normal HomoSapiens_BreastCancer_Normal --tumor HomoSapiens_BreastCancer_Tumor --pair BreastCancer --mode Intergenic
gammatext-scope decay-dfa     --pair BreastCancer --mode Intergenic
```

`dotnet test` runs three layers of regression tests: *manuscript anchors*
(the numbers printed in the article are asserted against `expected/`),
*pipeline golden tests* (`decay-dfa` always; `classify`, `density` and
`alpha-compare` on chromosome 21 when the sequence data is present) and chart
rendering smoke tests. Status of the lift against the legacy outputs
(chromosome 21, both cancers): `classify` reproduces the legacy window
classification exactly; `density` and `alpha-compare` were byte-identical
until the correction of 29 October 2026 described below.

Provenance of `expected/v1`: Intergenic and Gene outputs are the legacy runs
of 5 June / 15 September 2026; the `WholeGenomeMerged` and `WholeGenome`
tables were regenerated on 12 October 2026 with this code from the legacy
classifications (the legacy whole-genome run of 22 May predated the annotation
change of 27 May). `expected/v2` was produced entirely with this code on
15 October 2026 from the inputs built by `scripts/inputs/`.

Correction of 29 October 2026: the tumour count of a window longer than
200 kb (one that spans an N-run, a centromere or an acrocentric arm; about
1.5 % of all windows) was counted in 100-kb chunks, and a chunk whose
stop-codon density differed by more than 2x from the window mean was dropped
and the rest scaled up to the full window. An N-run has density 0, so its
chunks were always dropped and the real chunks were multiplied by 1.5, 2, 4 ...;
every |Delta| above about 100 in the `alpha-compare` tables of both datasets,
and every `Tumor_StopCodons = 0` row, came from that scaling and not from the
sequences. The rule now counts every chunk. The `StopCodonAlphaComparison_*`
tables of `expected/v1` and `expected/v2` were regenerated on that date; only
the `Tumor_StopCodons` and `Delta` columns of windows longer than 200 kb
changed, and the window classes, decay dots and `DecayGapDFA_*` summaries are
identical. The tables as they stood before are in the git history, and the
figures of the submitted manuscript show the artefact as off-scale peaks.

Runtime: `classify` is single-threaded per chromosome and its legacy window
scanner is quadratic in chromosome length (it re-tokenises the remainder of
the chromosome for every window), so chr21 takes ~1–2 min but chr1 ~50 min and
a genome ~5–6 h in one process. Run chromosomes as parallel processes instead
(`scripts/classify_parallel.sh`; one `--chromosomes` value each, ~2.5 GB RAM
per process) and a sample finishes in about an hour. The three experiments
take seconds per chromosome.

## Figures

`figures/make_figures.py` renders the article figures from a reference dataset
in `expected/` (matplotlib only):

```bash
python figures/make_figures.py --dataset v2      # -> figures/out/v2/
```

It produces Figure 1 (per-chromosome stop-codon difference; the breast-cancer
panel alone and both cancers side by side), Figures 2–3 (chromosome 1
DNA-scopes), Appendices 1–2 (plaited DNA-scopes of all chromosomes, npc-DNA
and pc-DNA, three A4 pages each) and `tables.md` with the DFA summaries. The
decay-dot rule follows `DecayGapDfa` (the dot counts per chromosome equal
the DFA summaries). The centromere band is the GRCh38 centromere
(`CentromereRegions.cs`) mapped into the window index through the normal
sample's gene mask with the same interval logic as `FilteredToRawMapper`;
this needs the sample data under `--data-root` (or `GAMMATEXT_DATA`) and
falls back to the legacy zero-stop-window heuristic without it.
`figures/out/v1/` and `figures/out/v2/` hold the rendered sets for the two
datasets; the v1 set reproduces the content of the submitted figures except
for the off-scale peaks, which were the counting artefact corrected on
29 October 2026 (see *Reproducing the article*).

## Input provenance

The sequences are pseudo-haploid consensus genomes built from public variant
data on GRCh38; the annotation is GENCODE v44 (Ensembl 110) lifted onto each
consensus. `scripts/inputs/` contains the scripts exactly as run and a README
with the construction rule.

| | HCC1395 / HCC1395BL (breast cancer) | COLO829 / COLO829BL (melanoma) |
|---|---|---|
| normal variants | `paftools.js call` (minimap2 2.30) of NCBI assembly **GCA_021234545.1** (HCC1395BL_v1.0) against GRCh38 | ONT open data, **Clair3 1.0.4** germline calls |
| tumour variants | SEQC2 high-confidence somatic **sSNV + sINDEL v1.2.1** (Fang et al. 2021) | ONT wf-somatic-variation, **ClairS 0.1.7**, PASS |
| tumour VCF | germline minus somatic positions, plus somatic sites | the normal's germline records (same order) minus those overlapping a somatic span, plus somatic |
| consensus | `bcftools consensus -H A` 1.13 per chromosome, both samples (ALT allele of every genotype) | same |
| annotation liftover | `minimap2 -x asm5` consensus vs GRCh38, `paftools.js liftover -q 0 -l 1000 -d 1.1`, per sample | same |
| consistency check | normal–tumour alignment breakpoints on chr21 = 21 (20 somatic indels) | 23 (22 somatic indels) |

The dataset-v1 inputs (submitted manuscript) differed from this in three ways
that are documented in `scripts/inputs/README.md`: the melanoma normal was
built without heterozygous alleles while the tumour carried them; the melanoma
tumour union ordered multi-allelic records differently from the normal, so
`bcftools consensus` applied different alleles at ~1,200 indel sites per
chromosome 21; and both breast-cancer samples used one annotation lifted, with
the alignment direction inverted, onto an earlier tumour FASTA in a coordinate
system up to ~100 kb off the analysed strings. The v2 construction removes all
three; the consistency checks above are part of the scripts.

## Differences from the legacy implementation

* All paths are configurable; nothing is hard-coded to the authors' drives.
* The window scanner tokenises incrementally instead of splitting the whole
  remainder of the chromosome for every window. It is equivalence-tested
  against the legacy scanner (`WindowScannerTests`) and reproduces the legacy
  classification files of every chromosome it has been compared with.
* The per-window PDF charts and `_DictionaryStats` dumps of the classification
  step are not produced (they are not consumed by any experiment).
* In the classification index files the lexicon size and the per-window
  elapsed time are now tab-separated (the legacy writer concatenated them).
  Columns 1–5, the only ones read downstream, are unchanged.
* Chart titles and legends are in English; the difference chart is named
  `Chart_<pair>_StopCodon_Difference.pdf`.
* `StopCodonStatistics` (mean negative/absolute delta) is not included; the
  article does not use it.
* Classification values in the tables keep the legacy spellings `alfa`,
  `zipf`, `ikke` so that outputs stay byte-comparable with `expected/`.
* Every output file is written with CRLF line endings and invariant-culture
  number formatting on every platform, so a run on Linux or on a machine with a
  non-English locale stays byte-identical with `expected/`.

## Licence

GPL-3.0-or-later. The DFA implementation is a translation of PhysioNet's
`dfa.c` (GPL-2.0-or-later); see `THIRD_PARTY_NOTICES.md` for the required
citations.
