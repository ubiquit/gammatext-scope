#!/usr/bin/env bash
# Melanoma v3 driver (Git Bash): WSL preparation -> pipeline -> verification. No backup
# needed (v1 backed up on 2026-10-15; v2 kept in place). Nothing existing is modified.
set -uo pipefail
export MSYS_NO_PATHCONV=1
V2="D:/Cancer sequences/v2"; LOG="$V2/run_v3.log"; STREAM="D:/Files/Stream"; HARV="D:/Files/Harvests"
REPO="D:/Source/gammatext-scope"; BIN="C:/tmp/gts-bin-v2"; EXE="$BIN/gammatext-scope.exe"; export EXE
N="HomoSapiens_Melanoma_Normal_v3"; T="HomoSapiens_Melanoma_Tumor_v3"
exec > >(tee -a "$LOG") 2>&1
stamp() { echo "[$(date +%H:%M:%S)] $*"; }
fail() { stamp "FAILED: $*"; exit 1; }

stamp "=== B. WSL: VCFs, consensus, annotation ==="
wsl.exe -e bash -lc "bash '/mnt/d/Cancer sequences/v2/prepare_v3.sh'" || fail "prepare_v3.sh"

stamp "=== C. pipeline ==="
[ -x "$EXE" ] || { rm -rf "$BIN"; cp -r "$REPO/src/GammaTextScope/bin/Release/net8.0" "$BIN"; }
for s in $N $T; do
  "$EXE" prepare-annotation --sample "$s" || fail "prepare-annotation $s"
  for m in Intergenic Gene WholeGenome; do
    stamp "classify $s $m"; bash "$REPO/scripts/classify_parallel.sh" "$s" "$m" 12 >/dev/null || fail "classify $s $m"
    n=$(ls "$HARV/$s/Experiment_${s}_${m}_7000/WindowWidth1"/*_indexLength.txt | wc -l); [ "$n" -eq 25 ] || fail "classify $s $m: $n/25"
  done
done
for m in Intergenic Gene WholeGenome; do
  "$EXE" density --normal $N --tumor $T --pair Melanoma_v3_Normal_vs_Tumor --mode $m > "$V2/density_Melanoma_v3_$m.log" 2>&1 || fail "density $m"
  "$EXE" alpha-compare --normal $N --tumor $T --pair Melanoma_v3 --mode $m > /dev/null 2>&1 || fail "alpha-compare $m"
  "$EXE" decay-dfa --pair Melanoma_v3 --mode $m > /dev/null 2>&1 || fail "decay-dfa $m"
done
"$EXE" alpha-compare --normal $N --tumor $T --pair Melanoma_v3 --mode WholeGenome --merged > /dev/null 2>&1 || fail "alpha-compare merged"
"$EXE" decay-dfa --pair Melanoma_v3 --mode WholeGenome --merged > /dev/null 2>&1 || fail "decay-dfa merged"

stamp "=== D. verification ==="
awk '/StopCodonDensity \(WholeGenome\): (21|1)\r?$/{c=$NF} /alignment:/{if(c!="")printf "chr%s %s\n", c, $0; c=""}' "$V2/density_Melanoma_v3_WholeGenome.log" | tr -d '\r'
python "$V2/verify_v2.py" v3 || fail "verify"
stamp "=== DFA summaries: v3 | v2 | v1 ==="
for m in Intergenic Gene WholeGenome WholeGenomeMerged; do
  echo "### Melanoma $m"; paste <(cut -f1,2,7 "$HARV/Tumor/Melanoma_v3/DecayGapDFA_$m/DecayGapDFA_Melanoma_v3_Summary.txt") <(cut -f2,7 "$HARV/Tumor/Melanoma_v2/DecayGapDFA_$m/DecayGapDFA_Melanoma_v2_Summary.txt" | sed 's/^/| /') <(cut -f2,7 "$REPO/expected/Tumor/Melanoma/DecayGapDFA_$m/DecayGapDFA_Melanoma_Summary.txt" | sed 's/^/| /') | column -t -s$'\t'
done
stamp "=== run_v3 finished ==="
