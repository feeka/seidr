#!/usr/bin/env bash
# Clone (if absent) and build MEGAHIT. Produces: megahit/build/megahit_core
# Requires: git, cmake, g++, zlib, make. Run under Linux/WSL.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

if [ ! -d megahit ]; then
  echo "[megahit] cloning..."
  git clone https://github.com/voutcn/megahit.git megahit
fi
echo "[megahit] configuring + building..."
cmake -S megahit -B megahit/build -DCMAKE_BUILD_TYPE=Release
cmake --build megahit/build -j"$(nproc)"
echo "[megahit] done -> $(./megahit/build/megahit_core dumpversion)"
