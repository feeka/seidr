#!/usr/bin/env bash
# Verify the §9.3 length gate REJECTS out-of-range ORFs. Plants 24 cassettes with
# lengths spanning <210, [210,1200], and >1200 bp; checks that exactly the in-range
# ones are recovered and the out-of-range ones are NOT.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
export PATH="$HERE/deps/bin:$PATH"
MH="$HERE/megahit/build/megahit_core"; B="$HERE/build"

# aa range [10,500] -> ORF lengths 3*aa+6 = 36..1506 bp (mix around the [210,1200] gate)
python3 sdbg_gen.py 3 "$B/reads.fa" "$B/truth.txt" "$B/vcr.fa" 24 10 500 2>/dev/null
printf "#lib\nse %s\n" "$B/reads.fa" > "$B/data.lib"
"$MH" buildlib "$B/data.lib" "$B/libpref" >/dev/null 2>&1
"$MH" read2sdbg --host_mem 2000000000 --read_lib_file "$B/libpref" \
      -m 1 -k 23 --num_cpu_threads 4 --o "$B/graph" >/dev/null 2>&1
"$B/vcr_array" "$B/graph" "$B/profile.txt" "$B/orfs.fa" 2>/dev/null

echo -n "planted ORF lengths: "
awk -F'\t' '/^orf/{print length($2)}' "$B/truth.txt" | sort -n | tr '\n' ' '; echo
python3 check_reject.py "$B/orfs.fa" "$B/truth.txt"
