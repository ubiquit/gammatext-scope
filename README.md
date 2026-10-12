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

`<mode>` is `Intergenic` (npc-DNA), `Gene` (pc-DNA) or `WholeGenome`
(Intergenic + Gene windows merged in genomic order).

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

`expected/` holds the pipeline outputs the manuscript was built from (tables,
DFA summaries and the per-window classifications of all four samples), so the
experiments can be re-run and diffed without the multi-hour classification step:

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
classification exactly; `density` and `alpha-compare` are byte-identical.

Provenance of `expected/`: Intergenic and Gene outputs are the legacy runs of
5 June / 15 September 2026. The legacy WholeGenome outputs (22 May 2026) were
produced from classification files that predate the 27 May annotation
harmonisation between tumour and normal; the WholeGenome tables shipped here
were regenerated on 12 October 2026 from the current classifications with this
code. One orphaned gap-series file from a March run was dropped.

Runtime (one core, 2026 desktop): `classify` ≈ 75 s for chromosome 21, i.e.
roughly 40 min per sample and mode genome-wide; the three experiments take
seconds per chromosome.

## Input provenance

The sequences are haplotype-1 consensus genomes built from public data on the
Ensembl GRCh38 primary assembly; the annotation is GENCODE v44 (Ensembl 110)
lifted onto each consensus. `scripts/prepare_inputs.sh` reconstructs them.

| | HCC1395 / HCC1395BL (breast cancer) | COLO829 / COLO829BL (melanoma) |
|---|---|---|
| normal variants | `paftools.js call` (minimap2 2.30) of NCBI assembly **GCA_021234545.1** (HCC1395BL_v1.0) against GRCh38 | ONT open data, **Clair3 1.0.4** germline calls |
| tumour variants | SEQC2 high-confidence somatic **sSNV + sINDEL v1.2.1** (Fang et al. 2021) | ONT wf-somatic-variation, **ClairS 0.1.7**, PASS |
| tumour VCF | germline minus somatic positions, plus somatic sites | germline minus somatic spans, plus somatic |
| consensus | `bcftools consensus` 1.13 per chromosome | `bcftools consensus` 1.13 per chromosome, with `-c` chain files |
| annotation liftover | `paftools.js liftover` on BED6 | CrossMap with the consensus chain files |

## Differences from the legacy implementation

* All paths are configurable; nothing is hard-coded to the authors' drives.
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

## Licence

GPL-3.0-or-later. The DFA implementation is a translation of PhysioNet's
`dfa.c` (GPL-2.0-or-later); see `THIRD_PARTY_NOTICES.md` for the required
citations.
