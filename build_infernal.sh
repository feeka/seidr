#!/usr/bin/env bash
# Build Infernal (gives cmsearch + cmemit) into ./deps  (no root required).
# Requires: curl, tar, gcc, make. Run under Linux/WSL.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
V=1.1.5

if [ ! -x "deps/bin/cmsearch" ]; then
  echo "[infernal] downloading + building $V..."
  mkdir -p deps/src && cd deps/src
  [ -f "infernal-$V.tar.gz" ] || curl -fsSL -o "infernal-$V.tar.gz" "http://eddylab.org/infernal/infernal-$V.tar.gz"
  rm -rf "infernal-$V" && tar xzf "infernal-$V.tar.gz"
  cd "infernal-$V"
  ./configure --prefix="$HERE/deps"
  make -j"$(nproc)"
  make install
fi
echo "[infernal] done -> $("$HERE/deps/bin/cmsearch" -h | grep INFERNAL | head -1)"
