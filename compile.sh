#!/usr/bin/env bash
# Compile the detector against MEGAHIT's SDBG sources -> build/vcr_fold.
# Requires build_megahit.sh to have been run.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
B="$HERE/build"; mkdir -p "$B"
MHSRC="$HERE/megahit/src"
SDBG_SRC="$MHSRC/sdbg/sdbg_meta.cpp $MHSRC/sdbg/sdbg_raw_content.cpp"

echo "[compile] vcr_fold"
g++ -O2 -std=c++17 -I "$MHSRC" "$HERE/vcr_fold.cpp" $SDBG_SRC -lpthread -o "$B/vcr_fold"
echo "[compile] done -> $B/vcr_fold"
