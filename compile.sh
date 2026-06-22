#!/usr/bin/env bash
# Compile our detector against MEGAHIT's SDBG sources. -> build/vcr_traverse
# (also build/vcr_anchor, build/sdbg_load_test). Requires build_megahit.sh run.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
B="$HERE/build"; mkdir -p "$B"
MHSRC="$HERE/megahit/src"
SDBG_SRC="$MHSRC/sdbg/sdbg_meta.cpp $MHSRC/sdbg/sdbg_raw_content.cpp"

for tool in vcr_traverse vcr_anchor sdbg_load_test; do
  echo "[compile] $tool"
  g++ -O2 -std=c++17 -I "$MHSRC" "$HERE/$tool.cpp" $SDBG_SRC -lpthread -o "$B/$tool"
done
echo "[compile] done -> $B/{vcr_traverse,vcr_anchor,sdbg_load_test}"
