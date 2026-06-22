#!/usr/bin/env python3
# Compare recovered ORFs (FASTA) against planted ORFs (truth.txt 'orf<i>\t<seq>').
# Prints a one-line verdict. Exit 0 iff every planted ORF was recovered byte-identical.
import sys
orfs_fa, truth_txt, label = sys.argv[1], sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "")
rec = set(l.strip() for l in open(orfs_fa) if not l.startswith(">") and l.strip())
planted = set()
for line in open(truth_txt):
    if line.startswith("orf"):
        planted.add(line.rstrip("\n").split("\t")[1])
match = len(planted & rec)
ok = (match == len(planted) and len(planted) > 0)
verdict = "OK" if ok else "FAIL"
print("%-18s planted=%d recovered=%d byte_identical=%d/%d extra=%d  %s"
      % (label, len(planted), len(rec), match, len(planted), len(rec - planted), verdict))
sys.exit(0 if ok else 1)
