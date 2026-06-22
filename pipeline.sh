#!/usr/bin/env bash
# Full pipeline, in order. Run under Linux/WSL from the repo root.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
bash build_megahit.sh      # 1. build MEGAHIT (clones it if absent)
bash build_infernal.sh     # 2. build Infernal (cmsearch + cmemit) into ./deps
bash generate_data.sh      # 3. make the VCR/profile/reads + build the SDBG
bash compile.sh            # 4. compile our detector against MEGAHIT's SDBG
bash run.sh                # 5. detect on the SDBG + validate with cmsearch
echo "=== pipeline complete; artifacts in ./build ==="
