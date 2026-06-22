#!/usr/bin/env python3
# Controlled synthetic super-integron for testing the SDBG anchor detection.
# Structure:  flank + VCR + (ORF + VCR) x N + flank
#   VCR identical across copies (so it collapses);
#   each ORF = ATG + non-stop codons + a RANDOMLY CHOSEN stop (TAA/TAG/TGA),
#   so the stop sits in the 3 nodes UPSTREAM of the (indeg>=2) VCR-start node.
import random, sys
seed   = int(sys.argv[1]) if len(sys.argv) > 1 else 1
out    = sys.argv[2] if len(sys.argv) > 2 else '/tmp/reads.fa'
truth  = sys.argv[3] if len(sys.argv) > 3 else '/tmp/truth.txt'
vcrfile= sys.argv[4] if len(sys.argv) > 4 else None   # real VCR (e.g. cmemit -c)
ncass  = int(sys.argv[5]) if len(sys.argv) > 5 else 6 # number of cassettes (ORFs)
aa_min = int(sys.argv[6]) if len(sys.argv) > 6 else 70   # ORF codon count range (len = 3*aa+6)
aa_max = int(sys.argv[7]) if len(sys.argv) > 7 else 390
random.seed(seed)
K = 23
flank = 100
readlen, stride = 250, 10
B = 'ACGT'; STOPS = ['TAA', 'TAG', 'TGA']
def rseq(n): return ''.join(random.choice(B) for _ in range(n))
def orf(aa):
    s = 'ATG'
    for _ in range(aa):
        c = rseq(3)
        while c in STOPS: c = rseq(3)
        s += c
    return s + random.choice(STOPS)            # varied stop
if vcrfile:                                     # use a REAL attC (so the CM profile applies)
    vcr = ''.join(l.strip() for l in open(vcrfile) if not l.startswith('>')).upper().replace('U','T')
else:
    vcr = rseq(80)
vcrlen = len(vcr)
arr = rseq(flank) + vcr
orfs = []
for _ in range(ncass):
    o = orf(random.randint(aa_min, aa_max))  # length 3*aa+6; default range ⊂ [210,1200]
    orfs.append(o)
    arr += o + vcr
arr += rseq(flank)
with open(out, 'w') as f:
    j = 0
    for p in range(0, len(arr) - readlen + 1, stride):
        f.write('>r%d\n%s\n' % (j, arr[p:p+readlen])); j += 1
    f.write('>r%d\n%s\n' % (j, arr[-readlen:])); j += 1
with open(truth, 'w') as f:
    f.write('vcr\t%s\n' % vcr)
    f.write('vcr_start_kmer\t%s\n' % vcr[:K])
    f.write('ncass\t%d\n' % ncass)
    for i, o in enumerate(orfs):                 # planted ORFs = ground truth for the beam
        f.write('orf%d\t%s\n' % (i, o))
sys.stderr.write('[gen] array_len=%d reads=%d vcr_start_23mer=%s\n' % (len(arr), j+1, vcr[:K]))
