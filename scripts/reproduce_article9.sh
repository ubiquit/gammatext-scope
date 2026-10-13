#!/usr/bin/env bash
# Regenerates every table used in the article from the per-window
# classifications in expected/ (or from freshly generated ones) and diffs the
# DFA summaries against expected/. Sequences must be present under $GAMMATEXT_DATA.
set -euo pipefail
EXE=${EXE:-src/GammaTextScope/bin/Release/net8.0/gammatext-scope}
OUT=${OUT:-./out}
export GAMMATEXT_OUT="$OUT"

# Stage the classifications (skip if you ran `classify` for all four samples yourself)
for s in HomoSapiens_BreastCancer_Normal HomoSapiens_BreastCancer_Tumor HomoSapiens_Melanoma_Normal HomoSapiens_Melanoma_Tumor; do
  for m in Intergenic Gene; do
    d="$OUT/$s/Experiment_${s}_${m}_7000/WindowWidth1"; mkdir -p "$d"
    cp -n expected/$s/Experiment_${s}_${m}_7000/WindowWidth1/*_window_positions_indexLength.txt "$d/" 2>/dev/null || true
  done
done

run_pair() {  # normal tumor density-pair alpha-pair
  for m in Intergenic Gene WholeGenome; do
    "$EXE" density --normal "$1" --tumor "$2" --pair "$3" --mode $m
  done
  for m in Intergenic Gene; do
    "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode $m
    "$EXE" decay-dfa     --pair "$4" --mode $m
  done
  # legacy whole-DNA variant: union of the Intergenic and Gene windows
  "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode WholeGenome --merged
  "$EXE" decay-dfa     --pair "$4" --mode WholeGenome --merged
  # annotation-free whole-DNA variant: needs `classify --mode WholeGenome` for both samples
  if "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode WholeGenome; then
    "$EXE" decay-dfa --pair "$4" --mode WholeGenome
  fi
}
run_pair HomoSapiens_BreastCancer_Normal HomoSapiens_BreastCancer_Tumor BreastCancer_Normal_vs_Tumor BreastCancer
run_pair HomoSapiens_Melanoma_Normal     HomoSapiens_Melanoma_Tumor     Melanoma_Normal_vs_Tumor     Melanoma

echo "== DFA summaries (Tables 2-3) vs expected/"
for p in BreastCancer Melanoma; do for m in Intergenic Gene WholeGenome; do
  f="Tumor/$p/DecayGapDFA_$m/DecayGapDFA_${p}_Summary.txt"
  if diff -q "$OUT/$f" "expected/$f" >/dev/null; then echo "identical  $f"; else echo "DIFFERS    $f"; fi
done; done
