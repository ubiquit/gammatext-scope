#!/usr/bin/env bash
set -euo pipefail
# Usage: ./make_canonical_json.sh <FAI> <OUT.json> "<NAME>" <ACC> [TAXID]

FAI="${1:?FAI index required}"
OUT="${2:?Output JSON path required}"
NAME="${3:?Assembly display name required}"
ACC="${4:?Assembly ACC code required}"
TAXID="${5:-9606}"

awk -v name="$NAME" -v acc="$ACC" -v taxid="$TAXID" '
  # Load lengths, normalize chr prefixes and MT/M
  FNR==NR{
    len[$1]=$2
    tmp=$1; gsub(/^chr/,"",tmp); len[tmp]=$2
    if ($1=="MT") len["M"]=$2
    if ($1=="M")  len["MT"]=$2
    next
  }
  END{
    n = split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT", order, " ")
    print "{"
    printf("  \"name\": \"%s\",\n", name)
    printf("  \"acc\": \"%s\",\n", acc)
    printf("  \"taxid\": %d,\n", taxid)
    printf("  \"type\": \"Primary-Unit\",\n")
    printf("  \"seq\": [\n")
    for (i=1; i<=n; i++){
      k = order[i]
      L = (k in len) ? len[k] : "null"
      printf("    {\"name\":\"%s\",\"chrom\":1,\"topl\":1,\"length\":%s,\"acc\":\"%s\",\"gi\":null}", k, L, k)
      if (i < n) printf(",")
      printf("\n")
    }
    print "  ]"
    print "}"
  }
' "$FAI" /dev/null > "$OUT"

echo "Wrote $OUT"
