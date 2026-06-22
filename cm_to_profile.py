#!/usr/bin/env python3
# Parse an Infernal .cm into a LINEAR per-position nucleotide profile (option a):
#   MATL/MATR match states -> the 4 singlet log-odds (A C G U) directly;
#   MATP pair states        -> marginalize the 16 pair log-odds to a singlet
#                              per column.
# Emissions are log2-odds vs a uniform null (NULL line = 0.000). Marginal of a
# pair at the LEFT column for base x:  log2( 0.25 * sum_y 2^pairLO[4x+y] ).
# Output: one line per consensus column:  pos  loA  loC  loG  loT   (U==T).
import sys, math
cm  = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else '/tmp/profile.txt'
lines = open(cm).read().splitlines()

clen = None
for l in lines:
    if l.startswith('CLEN'): clen = int(l.split()[1])
cur = lmap = rmap = None
cols = []                             # (MAP_coord, [loA,loC,loG,loT])
def nums(s):
    o = []
    for t in s.split():
        try: o.append(float(t))
        except ValueError: pass
    return o
for l in lines:
    s = l.strip()
    if not s: continue
    if s.startswith('['):             # node line: [ TYPE idx ] lmap rmap lchar rchar
        before, after = s.split(']', 1)
        cur = before.replace('[', '').split()[0]
        a = after.split()
        lmap = int(a[0]) if a[0] != '-' else None   # MAP/alignment coord, NOT 1..CLEN
        rmap = int(a[1]) if a[1] != '-' else None
        continue
    st = s.split()[0]
    if cur == 'MATL' and st == 'ML' and lmap is not None:
        cols.append((lmap, nums(s)[-4:]))
    elif cur == 'MATR' and st == 'MR' and rmap is not None:
        cols.append((rmap, nums(s)[-4:]))
    elif cur == 'MATP' and st == 'MP':
        em = nums(s)[-16:]            # AA AC AG AU CA ... UU
        L = [math.log2(0.25 * sum(2**em[4*x+y] for y in range(4))) for x in range(4)]
        R = [math.log2(0.25 * sum(2**em[4*x+y] for x in range(4))) for y in range(4)]
        if lmap is not None: cols.append((lmap, L))
        if rmap is not None: cols.append((rmap, R))

cols.sort(key=lambda t: t[0])         # linear consensus order = sorted by MAP coord
with open(out, 'w') as f:
    for i, (mp, v) in enumerate(cols, 1):
        f.write('%d\t%.4f\t%.4f\t%.4f\t%.4f\n' % (i, v[0], v[1], v[2], v[3]))

cons = ''.join('ACGT'[max(range(4), key=lambda k: v[k])] for mp, v in cols)
sys.stderr.write('parsed %d match columns (CLEN=%d) -> %s\n' % (len(cols), clen, out))
print(cons)                           # per-column argmax consensus
