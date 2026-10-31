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
unmasked chromosome); `Exon` is accepted but experimental and untested. For `alpha-compare` and `decay-dfa`, `--mode WholeGenome`
requires whole-genome window classifications for both samples; `--merged`
selects the legacy behaviour instead — the union of the Intergenic and Gene
windows in genomic order (outputs suffixed `WholeGenomeMerged`), which depends
on the gene annotation.

## Build and run

Requires the .NET SDK 8.0 or later.

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

Correction of 31 October 2026: the offset map (`AlignmentOffsetMap`) re-anchored
with its wide search only every 10 Mb, so after a long N-run (a centromere),
where the anchors match trivially on N, it could come out with an offset wrong
by more than its ±200 bp search radius and stay lost until the next 10-Mb
point. A base-by-base check of every window found this in two stretches of
the breast-cancer pair (chr1 windows 645–669 and chr9 windows 292–346 in
npc-DNA; 80 of 16,438 windows), where the tumour window therefore covered a
different piece of DNA and the Delta values were noise; all other windows of
both pairs and both modes covered the same DNA (2–4 mismatched bases per
window). The map now re-anchors wherever the fine scan has failed 20 steps
(10 kb) in a row. The `alpha-compare` and `density` tables of both datasets
were regenerated; only windows and 100-kb bins inside such stretches changed
(dataset v2: 81 npc windows of the breast-cancer pair and the chr1/chr9 bars of
its Figure 1 by 2 and 1; nothing in the melanoma pair). The v1 melanoma pair,
whose inputs differ at every heterozygous site, was misaligned in 19 % of its
windows before and is still misaligned in 15 % after the correction; its v1
tables are kept for the record only.

Correction of 1 November 2026: the `P` and `Q` columns of the
`*_100k_diff.tsv` and `*_100k_significant.tsv` density tables used the error
function in place of the normal CDF, so `P` was 2·erfc(|z|) instead of
erfc(|z|/√2) and too small by a factor growing with |z|. The formula was
corrected and those two files regenerated for both datasets; counts, rates
and z are unchanged, and the article uses only the counts.

Runtime: `classify` is single-threaded per chromosome; with the incremental
window scanner chr21 takes a few seconds and chr1 about 20 s, a whole sample a
few minutes in one process, less with `scripts/classify_parallel.sh` (one
process per chromosome, ~2.5 GB RAM each). `density` and `alpha-compare` load
one chromosome of both samples at a time (~1.5 GB RAM) and take a few seconds
per chromosome; `decay-dfa` runs in seconds.

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
| normal variants | `paftools.js call` (minimap2 2.24, paftools.js 2.30) of NCBI assembly **GCA_021234545.1** (HCC1395BL_v1.0) against GRCh38 | ONT open data, **Clair3 1.0.4** germline calls |
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

## Legacy conventions pinned by `expected/`

The following behaviours of the lifted code are deliberate reproductions of
the legacy platform. They are pinned by the byte-identity tests and must not
be "fixed" without regenerating `expected/`:

* **Gene mask offset.** The protein-coding mask is applied in an index space
  shifted by two: a GFF gene `[start, end]` masks the 0-based positions
  `start−2 … end−2` (`WindowHelper`, `FilteredToRawMapper`). The first
  (lowest-start) protein-coding gene of every chromosome is skipped; in Gene
  mode its bases are never emitted, in Intergenic mode they count as npc-DNA.
  A gene nested inside its predecessor leaves one base misassigned after the
  outer gene ends. `PipelineTests.GeneMaskKeepsTheLegacyOffsetAndSkipsTheFirstGene`
  guards this.
* **Windows.** A window is the shortest prefix of the stop-codon tokenisation
  with more than 7,000 tokens that ends in a non-empty token (so a window holds
  7,001 stop codons, 7,002 when a token is empty); the last window of a
  chromosome whose text ends with a stop codon is reported 3 bases too long.
* **Classification.** The wizard-hat histogram bins are built by repeated
  `+= 0.1` and compared with the literals 0.6, 0.7 and 1.0; this is
  deterministic and `ClassificationTests` guards it. Ties for the largest peak
  go to the lowest α bin. Class names are the legacy `alfa`, `zipf`, `ikke`.
* **Window pairing.** `alpha-compare` pairs the i-th window of the normal with
  the i-th window of the tumour: `StartPos`/`Length`/`Normal_StopCodons` are
  the normal's window, `Tumor_Class` is the tumour's own i-th window (which
  has a different extent), and `Tumor_StopCodons` is counted on the tumour
  range that the offset map assigns to the normal window. Surplus windows of
  the longer sample are dropped.
* **Decay dots and DFA.** A decay dot marks a 10-window region `[c−5, c+5)`
  in which one sample has at least 3 more qualified windows than the other;
  in the genome-wide series the region slides across chromosome boundaries;
  ties between alfa and Zipf dominance count as Zipf. DFA α is the global
  log–log slope over box sizes 4 … N/4 (no crossover selection) and is
  reported for 30 or more gaps; `N/A (DFA error)` marks series of 15–19 gaps,
  `N/A (insufficient data)` fewer than 15.
* **Density.** `density` counts stop codons on 100-kb bins of the normal
  mapped through the offset map; "trials" count triplet starts including N,
  and the genome-summary `BinEnd` is a cumulative length, not a coordinate.
* **Alignment.** The offset map treats N as a wildcard, so it is blind inside
  assembly gaps and re-anchors after them (see the correction of 31 October
  above). It assumes the two sequences share their germline; for a pair with
  a difference every kilobase (the v1 melanoma inputs) it loses the sequence
  in a sizeable share of windows.
* **Input contract.** One file per chromosome, upper-case `ACGTN`, no FASTA
  header, named `chrN.<anything>` (the chromosome is the last `_`-separated
  field before the first `.`); other files in a sample folder are ignored with
  a warning.

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

## How to cite

See `CITATION.cff`. Please cite the article above and, for the software, the
Zenodo record of the release you used (DOI to be added with the first
release). The DFA implementation requires the PhysioNet citations listed in
`THIRD_PARTY_NOTICES.md`.

## Licence

GPL-3.0-or-later. The DFA implementation is a translation of PhysioNet's
`dfa.c` (GPL-2.0-or-later); see `THIRD_PARTY_NOTICES.md` for the required
citations.
