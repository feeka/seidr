#!/usr/bin/env bash
# Stress-test ORF recovery against planted ground truth across several seeds and
# cassette counts. Each dataset: regenerate reads (distinct ORFs, identical VCR) ->
# rebuild the SDBG -> run vcr_array -> compare recovered vs planted byte-for-byte.
# Requires generate_data.sh + compile.sh already run once (needs build/profile.txt,
# build/vcr.fa, build/vcr_array, and the built MEGAHIT/Infernal deps).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
export PATH="$HERE/deps/bin:$PATH"
MH="$HERE/megahit/build/megahit_core"; B="$HERE/build"

# datasets: "seed ncass" pairs
DATASETS=("1 6" "2 6" "3 6" "4 6" "5 6" "7 12" "9 20" "11 50")

fails=0
for ds in "${DATASETS[@]}"; do
  set -- $ds; seed=$1; nc=$2
  python3 sdbg_gen.py "$seed" "$B/reads.fa" "$B/truth.txt" "$B/vcr.fa" "$nc" 2>/dev/null
  printf "#lib\nse %s\n" "$B/reads.fa" > "$B/data.lib"
  "$MH" buildlib "$B/data.lib" "$B/libpref" >/dev/null 2>&1
  "$MH" read2sdbg --host_mem 2000000000 --read_lib_file "$B/libpref" \
        -m 1 -k 23 --num_cpu_threads 4 --o "$B/graph" >/dev/null 2>&1
  "$B/vcr_array" "$B/graph" "$B/profile.txt" "$B/orfs.fa" 2>/dev/null
  python3 check_orf.py "$B/orfs.fa" "$B/truth.txt" "seed=$seed ncass=$nc" || fails=$((fails+1))
done

echo "----"
if [ "$fails" -eq 0 ]; then echo "ALL DATASETS PASSED"; else echo "$fails dataset(s) FAILED"; fi
exit "$fails"
