#!/usr/bin/env bash
# Generate the test data and build the SDBG, all into ./build:
#   - a real attC VCR (cmemit -c)         -> build/vcr.fa
#   - the linear CM profile               -> build/profile.txt
#   - synthetic reads (VCR + ORFs)        -> build/reads.fa
#   - the MEGAHIT succinct de Bruijn graph-> build/graph.*
# Requires: build_infernal.sh and build_megahit.sh already run.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
export PATH="$HERE/deps/bin:$PATH"
MH="$HERE/megahit/build/megahit_core"
CM="$HERE/models/Vibrionales.cm"
B="$HERE/build"; mkdir -p "$B"

echo "[data] emitting a real attC VCR (cmemit -c)"
cmemit -c "$CM" | grep -v '^>' | tr -d '\n' | tr 'Uu' 'Tt' | tr 'a-z' 'A-Z' > "$B/vcr.seq"
printf ">vcr\n%s\n" "$(cat "$B/vcr.seq")" > "$B/vcr.fa"

echo "[data] CM -> linear per-position profile"
python3 cm_to_profile.py "$CM" "$B/profile.txt" >/dev/null

echo "[data] building synthetic reads (flank + VCR + (ORF+VCR)xN + flank)"
python3 sdbg_gen.py 1 "$B/reads.fa" "$B/truth.txt" "$B/vcr.fa"

echo "[data] MEGAHIT: buildlib + read2sdbg (k=23, -m 1) -> build/graph.*"
printf "#lib\nse %s\n" "$B/reads.fa" > "$B/data.lib"
"$MH" buildlib "$B/data.lib" "$B/libpref" >/dev/null 2>&1
"$MH" read2sdbg --host_mem 2000000000 --read_lib_file "$B/libpref" \
      -m 1 -k 23 --num_cpu_threads 4 --o "$B/graph" >/dev/null 2>&1
echo "[data] done. SDBG files:"; ls "$B"/graph.sdbg*
