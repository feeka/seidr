#!/usr/bin/env bash
# Run detection on the SDBG and validate:
#   build/vcr_traverse build/graph build/profile.txt -> build/cands.fa
#   cmsearch validates candidates; cmp checks recovery against the planted VCR.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
export PATH="$HERE/deps/bin:$PATH"
B="$HERE/build"
CM="$HERE/models/Vibrionales.cm"

echo "[run] traversing the SDBG (anchors -> greedy profile-guided walk)"
"$B/vcr_traverse" "$B/graph" "$B/profile.txt" "$B/cands.fa"
echo "[run] candidates:"; cat "$B/cands.fa"

echo "[run] cmsearch validation"
cmsearch --tblout "$B/cands.tbl" --noali -E 1 "$CM" "$B/cands.fa" >/dev/null 2>&1
grep -v '^#' "$B/cands.tbl" | awk '{print "  hit "$1" score="$15" E="$16}' || echo "  (no hits)"

echo "[run] recovery check vs planted VCR (cmp)"
grep -v '^>' "$B/cands.fa" | tr -d '\n' > "$B/cand.seq"
if cmp -s "$B/vcr.seq" "$B/cand.seq"; then echo "  BYTE-IDENTICAL to planted VCR"
else echo "  differs:"; cmp "$B/vcr.seq" "$B/cand.seq" || true; fi

echo "[run] ORF crossing (vcr_array: VCR -> START -> bounded beam) -> build/orfs.fa"
"$B/vcr_array" "$B/graph" "$B/profile.txt" "$B/orfs.fa"

echo "[run] ORF recovery check vs planted ground truth (set comparison)"
grep '^orf' "$B/truth.txt" | cut -f2 | sort -u > "$B/planted_orfs.txt"
grep -v '^>' "$B/orfs.fa" | sort -u                > "$B/recovered_orfs.txt"
np=$(wc -l < "$B/planted_orfs.txt"); nr=$(wc -l < "$B/recovered_orfs.txt")
match=$(comm -12 "$B/planted_orfs.txt" "$B/recovered_orfs.txt" | wc -l)
echo "  planted=$np  recovered=$nr  exact_matches=$match"
if [ "$np" -gt 0 ] && [ "$match" -eq "$np" ]; then
  echo "  ALL $np planted ORFs recovered BYTE-IDENTICAL"
else
  echo "  MISMATCH — planted ORFs NOT all recovered. Missing:"
  comm -23 "$B/planted_orfs.txt" "$B/recovered_orfs.txt" | sed 's/^/    /' | head
fi
