#!/usr/bin/env bash
# Regenerates every table of a reference dataset from the per-window classifications
# in expected/<dataset> (or from freshly generated ones) and diffs the DFA summaries.
# Sequences must be present under $GAMMATEXT_DATA.
#   DATASET=v2 scripts/reproduce_article9.sh   (default; v1 = submitted manuscript)
set -euo pipefail
DATASET=${DATASET:-v2}
EXE=${EXE:-src/GammaTextScope/bin/Release/net8.0/gammatext-scope}
OUT=${OUT:-./out/$DATASET}
export GAMMATEXT_OUT="$OUT"
EXP="expected/$DATASET"

if [ "$DATASET" = "v1" ]; then
  BC_N=HomoSapiens_BreastCancer_Normal; BC_T=HomoSapiens_BreastCancer_Tumor; BC_D=BreastCancer_Normal_vs_Tumor; BC_A=BreastCancer
  ML_N=HomoSapiens_Melanoma_Normal;     ML_T=HomoSapiens_Melanoma_Tumor;     ML_D=Melanoma_Normal_vs_Tumor;     ML_A=Melanoma
  MODES="Intergenic Gene"
else
  BC_N=HomoSapiens_BreastCancer_Normal_v2; BC_T=HomoSapiens_BreastCancer_Tumor_v2; BC_D=BreastCancer_v2_Normal_vs_Tumor; BC_A=BreastCancer_v2
  ML_N=HomoSapiens_Melanoma_Normal_v3;     ML_T=HomoSapiens_Melanoma_Tumor_v3;     ML_D=Melanoma_v3_Normal_vs_Tumor;     ML_A=Melanoma_v3
  MODES="Intergenic Gene WholeGenome"
fi

# Stage the classifications (skip if you ran `classify` for all four samples yourself)
for s in $BC_N $BC_T $ML_N $ML_T; do
  for m in $MODES; do
    d="$OUT/$s/Experiment_${s}_${m}_7000/WindowWidth1"; mkdir -p "$d"
    cp -n "$EXP/$s/Experiment_${s}_${m}_7000/WindowWidth1/"*_window_positions_indexLength.txt "$d/" 2>/dev/null || true
  done
done

run_pair() {  # normal tumor density-pair alpha-pair
  for m in Intergenic Gene WholeGenome; do "$EXE" density --normal "$1" --tumor "$2" --pair "$3" --mode $m; done
  for m in $MODES; do
    "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode $m
    "$EXE" decay-dfa     --pair "$4" --mode $m
  done
  "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode WholeGenome --merged
  "$EXE" decay-dfa     --pair "$4" --mode WholeGenome --merged
}
run_pair $BC_N $BC_T $BC_D $BC_A
run_pair $ML_N $ML_T $ML_D $ML_A

echo "== DFA summaries vs $EXP"
for p in $BC_A $ML_A; do for m in $MODES WholeGenomeMerged; do
  f="Tumor/$p/DecayGapDFA_$m/DecayGapDFA_${p}_Summary.txt"
  if diff -q "$OUT/$f" "$EXP/$f" >/dev/null; then echo "identical  $f"; else echo "DIFFERS    $f"; fi
done; done
