#!/usr/bin/env bash
# Classify every chromosome of a sample as its own process, N at a time.
# usage: scripts/classify_parallel.sh <sample> <Intergenic|Gene|WholeGenome> [parallel=8] [chromosomes...]
# Chromosomes already present in the output folder are skipped, so the script can resume.
set -euo pipefail
SAMPLE=${1:?sample}; MODE=${2:?mode}; PAR=${3:-8}; shift 3 2>/dev/null || shift $#
EXE=${EXE:-$(dirname "$0")/../src/GammaTextScope/bin/Release/net8.0/gammatext-scope}
OUT=${GAMMATEXT_OUT:-D:/Files/Harvests}
DATA=${GAMMATEXT_DATA:-D:/Files/Stream}
DEST="$OUT/$SAMPLE/Experiment_${SAMPLE}_${MODE}_7000/WindowWidth1"

if [ $# -gt 0 ]; then CHRS="$*"; else
  CHRS=$(ls "$DATA/$SAMPLE" | sed -E 's/^chr([^.]+)\..*/\1/' | tr '\n' ' ')
fi

running=0
for c in $CHRS; do
  if [ -f "$DEST/chr${c}WW1_window_positions_indexLength.txt" ]; then echo "skip chr$c (done)"; continue; fi
  ( t0=$(date +%s); "$EXE" classify --sample "$SAMPLE" --mode "$MODE" --chromosomes "$c" >/dev/null 2>&1 \
      && echo "chr$c done in $(( $(date +%s) - t0 ))s" || echo "chr$c FAILED" ) &
  running=$((running + 1))
  if [ "$running" -ge "$PAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "$SAMPLE $MODE: $(ls "$DEST"/*_indexLength.txt | wc -l) chromosome files in $DEST"
