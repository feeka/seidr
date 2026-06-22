#!/usr/bin/env python3
# Parse an Infernal .cm into a PAIR MODEL for idea-B (structure-aware) search.
# Unlike cm_to_profile.py (which FLATTENS pairs to singlets), this keeps the
# base-PAIR emissions and the consensus secondary structure.
#
# Output: one line per consensus column, 5'->3' order (renumbered 1..N):
#   col  role  partner  v0 v1 v2 ...
#     SL / SR : single (MATL/MATR). partner=-1. 4 singlet log-odds (A C G T).
#     P5      : 5' half of a pair.  partner=<3' col>. 4 MARGINAL-LEFT log-odds.
#     P3      : 3' half of a pair.  partner=<5' col>. 16 CONDITIONAL log-odds,
#               cond[x*4+y] = pairLO(x,y) - marginalLeft(x).
#   So score(5' base x via P5 marginal) + score(3' base y via P3 cond[x][y])
#   = full pair log-odds pairLO(x,y), with no double counting.
# Emissions are log2-odds vs a uniform null (as in the .cm MP/ML/MR lines).
import sys, math
cm  = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else '/tmp/pairmodel.txt'
lines = open(cm).read().splitlines()

clen = None
for l in lines:
    if l.startswith('CLEN'): clen = int(l.split()[1])

def nums(s):
    o = []
    for t in s.split():
        try: o.append(float(t))
        except ValueError: pass
    return o

cur = lmap = rmap = None
singles = {}    # map_coord -> (role 'SL'/'SR', [4 log-odds])
pairs   = []    # (lmap, rmap, [16 log-odds])  order AA AC AG AU CA ... UU
for l in lines:
    s = l.strip()
    if not s: continue
    if s.startswith('['):                       # node line: [ TYPE idx ] lmap rmap ...
        before, after = s.split(']', 1)
        cur = before.replace('[', '').split()[0]
        a = after.split()
        lmap = int(a[0]) if a[0] != '-' else None
        rmap = int(a[1]) if a[1] != '-' else None
        continue
    st = s.split()[0]
    if   cur == 'MATL' and st == 'ML' and lmap is not None: singles[lmap] = ('SL', nums(s)[-4:])
    elif cur == 'MATR' and st == 'MR' and rmap is not None: singles[rmap] = ('SR', nums(s)[-4:])
    elif cur == 'MATP' and st == 'MP':                       pairs.append((lmap, rmap, nums(s)[-16:]))

flat = (len(sys.argv) > 3 and sys.argv[3] == 'flat')   # 'flat' = singlet-only control (no pairs)

col = {}    # map_coord -> (role, partner_mapcoord, [vals], [singlet4])
for lmap, rmap, em in pairs:
    L = [math.log2(0.25 * sum(2**em[4*x+y] for y in range(4))) for x in range(4)]   # marginal-left
    R = [math.log2(0.25 * sum(2**em[4*x+y] for x in range(4))) for y in range(4)]   # marginal-right
    cond = [em[4*x+y] - L[x] for x in range(4) for y in range(4)]                    # conditional 3'
    if lmap is not None: col[lmap] = ('P5', rmap, L, L)
    if rmap is not None: col[rmap] = ('P3', lmap, cond, R)
for mp, (role, vals) in singles.items():
    col[mp] = (role, -1, vals, vals)

order = sorted(col.keys())                       # 5'->3' consensus order
idx = {mp: i + 1 for i, mp in enumerate(order)}  # renumber 1..N
with open(out, 'w') as f:
    for mp in order:
        role, partner, vals, sgl = col[mp]
        if flat:                                 # all columns -> singlet (marginal), no pairing
            f.write('%d\tSL\t-1\t%s\n' % (idx[mp], '\t'.join('%.4f' % v for v in sgl)))
        else:
            p = idx[partner] if (partner != -1 and partner in idx) else -1
            f.write('%d\t%s\t%d\t%s\n' % (idx[mp], role, p, '\t'.join('%.4f' % v for v in vals)))

sys.stderr.write('parsed %d columns (%d pairs, %d singles, CLEN=%s) -> %s\n'
                 % (len(order), len(pairs), len(singles), clen, out))
