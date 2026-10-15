#!/usr/bin/env bash
# Driver for the v2 regeneration (run from Git Bash on Windows):
#   A. backup of the v1 inputs and results      -> D:\Files\Backups\2026-10-15_article9_v1
#   B. WSL: consensus genomes + annotation (prepare_v2.sh) -> D:\Files\Stream\*_v2*
#   C. gammatext-scope: prepare-annotation, classify (3 modes), density, alpha-compare, decay-dfa
#   D. verification
# Everything is written to new *_v2 names; nothing existing is modified or deleted.
set -uo pipefail
export MSYS_NO_PATHCONV=1   # keep /E /MT:16 etc. as robocopy switches under Git Bash
V2="D:/Cancer sequences/v2"; LOG="$V2/run_v2.log"
STREAM="D:/Files/Stream"; HARV="D:/Files/Harvests"; BK="D:/Files/Backups/2026-10-15_article9_v1"
REPO="D:/Source/gammatext-scope"; BIN="C:/tmp/gts-bin-v2"
SAMPLES="HomoSapiens_BreastCancer_Normal HomoSapiens_BreastCancer_Tumor HomoSapiens_Melanoma_Normal HomoSapiens_Melanoma_Tumor"
exec > >(tee -a "$LOG") 2>&1
stamp() { echo "[$(date +%H:%M:%S)] $*"; }
fail() { stamp "FAILED: $*"; exit 1; }
w() { echo "$1" | sed 's|/|\\|g'; }   # D:/x -> D:\x for robocopy

if [ -f "$BK/.complete" ]; then stamp "=== A. backup already complete, skipping ==="; else
stamp "=== A. backup ==="
mkdir -p "$BK/Stream" "$BK/Harvests/Tumor"
for s in $SAMPLES; do
  for suf in "" _GeneInfo _AccessionProcessing _CommonWords; do
    [ -d "$STREAM/$s$suf" ] || continue
    robocopy "$(w "$STREAM/$s$suf")" "$(w "$BK/Stream/$s$suf")" /E /MT:16 /NFL /NDL /NJH /NJS /NP >/dev/null; rc=$?; [ $rc -ge 8 ] && fail "robocopy $s$suf rc=$rc"
  done
  for m in Intergenic Gene WholeGenome; do
    src="$HARV/$s/Experiment_${s}_${m}_7000/WindowWidth1"; [ -d "$src" ] || continue
    robocopy "$(w "$src")" "$(w "$BK/Harvests/$s/Experiment_${s}_${m}_7000/WindowWidth1")" "*_window_positions_indexLength.txt" /MT:8 /NFL /NDL /NJH /NJS /NP >/dev/null; rc=$?; [ $rc -ge 8 ] && fail "robocopy harvest $s $m rc=$rc"
  done
done
for p in BreastCancer BreastCancer_Normal_vs_Tumor Melanoma Melanoma_Normal_vs_Tumor; do
  robocopy "$(w "$HARV/Tumor/$p")" "$(w "$BK/Harvests/Tumor/$p")" /E /MT:16 /NFL /NDL /NJH /NJS /NP >/dev/null; rc=$?; [ $rc -ge 8 ] && fail "robocopy Tumor/$p rc=$rc"
done
touch "$BK/.complete"; stamp "backup done: $(du -sh "$BK" | cut -f1)"
fi

stamp "=== B. WSL: consensus + annotation ==="
wsl.exe -e bash -lc "bash '/mnt/d/Cancer sequences/v2/prepare_v2.sh'" || fail "prepare_v2.sh"

stamp "=== C. pipeline ==="
rm -rf "$BIN"; cp -r "$REPO/src/GammaTextScope/bin/Release/net8.0" "$BIN"; EXE="$BIN/gammatext-scope.exe"
export EXE
for s in $SAMPLES; do
  "$EXE" prepare-annotation --sample "${s}_v2" || fail "prepare-annotation ${s}_v2"
  stamp "${s}_v2 GeneInfo: $(wc -l < "$STREAM/${s}_v2_GeneInfo/${s}_v2_Annotation_Converted.txt") lines"
done
for s in $SAMPLES; do for m in Intergenic Gene WholeGenome; do
  stamp "classify ${s}_v2 $m"
  bash "$REPO/scripts/classify_parallel.sh" "${s}_v2" "$m" 12 >/dev/null || fail "classify ${s}_v2 $m"
  n=$(ls "$HARV/${s}_v2/Experiment_${s}_v2_${m}_7000/WindowWidth1"/*_indexLength.txt | wc -l); [ "$n" -eq 25 ] || fail "classify ${s}_v2 $m produced $n/25 files"
done; done
run_pair() {  # normal tumor density-pair alpha-pair
  for m in Intergenic Gene WholeGenome; do
    "$EXE" density --normal "$1" --tumor "$2" --pair "$3" --mode $m > "$V2/density_$3_$m.log" 2>&1 || fail "density $3 $m"
  done
  for m in Intergenic Gene WholeGenome; do
    "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode $m > "$V2/alpha_$4_$m.log" 2>&1 || fail "alpha-compare $4 $m"
    "$EXE" decay-dfa --pair "$4" --mode $m > /dev/null 2>&1 || fail "decay-dfa $4 $m"
  done
  "$EXE" alpha-compare --normal "$1" --tumor "$2" --pair "$4" --mode WholeGenome --merged > /dev/null 2>&1 || fail "alpha-compare $4 merged"
  "$EXE" decay-dfa --pair "$4" --mode WholeGenome --merged > /dev/null 2>&1 || fail "decay-dfa $4 merged"
  stamp "experiments done for $4"
}
run_pair HomoSapiens_BreastCancer_Normal_v2 HomoSapiens_BreastCancer_Tumor_v2 BreastCancer_v2_Normal_vs_Tumor BreastCancer_v2
run_pair HomoSapiens_Melanoma_Normal_v2     HomoSapiens_Melanoma_Tumor_v2     Melanoma_v2_Normal_vs_Tumor     Melanoma_v2

stamp "=== D. verification ==="
for p in BreastCancer_v2 Melanoma_v2; do
  echo "$p WholeGenome alignment breakpoints (chr21, chr1):"; grep -E "StopCodonDensity \(WholeGenome\): (21|1)$" -A1 "$V2/density_${p}_Normal_vs_Tumor_WholeGenome.log" | grep alignment
done
python "$V2/verify_v2.py" || fail "verify_v2.py"
stamp "=== DFA summaries: v2 vs v1 ==="
for p in BreastCancer Melanoma; do for m in Intergenic Gene WholeGenome WholeGenomeMerged; do
  f2="$HARV/Tumor/${p}_v2/DecayGapDFA_$m/DecayGapDFA_${p}_v2_Summary.txt"; f1="$REPO/expected/Tumor/$p/DecayGapDFA_$m/DecayGapDFA_${p}_Summary.txt"
  echo "### $p $m (v2 | v1)"; paste <(cut -f1,2,7 "$f2") <(cut -f2,7 "$f1" | sed 's/^/| /') | column -t -s$'\t'
done; done
stamp "=== run_v2 finished ==="
