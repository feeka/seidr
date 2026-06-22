#!/usr/bin/env python3
# Verify the §9.3 length gate: planted ORFs with length in [210,1200] must be
# recovered; those outside must NOT be. Exit 0 iff the gate behaves exactly.
import sys
orfs_fa, truth_txt = sys.argv[1], sys.argv[2]
LMIN, LMAX = 210, 1200
rec = set(l.strip() for l in open(orfs_fa) if not l.startswith(">") and l.strip())
planted = []
for line in open(truth_txt):
    if line.startswith("orf"):
        planted.append(line.rstrip("\n").split("\t")[1])
in_range  = [s for s in planted if LMIN <= len(s) <= LMAX]
out_range = [s for s in planted if not (LMIN <= len(s) <= LMAX)]
in_ok  = sum(1 for s in in_range  if s in rec)
out_bad= sum(1 for s in out_range if s in rec)   # out-of-range that leaked through (should be 0)
ok = (in_ok == len(in_range)) and (out_bad == 0)
print("in-range planted=%d recovered=%d | out-of-range planted=%d wrongly_recovered=%d  %s"
      % (len(in_range), in_ok, len(out_range), out_bad, "OK" if ok else "FAIL"))
sys.exit(0 if ok else 1)
