"""Does GRCh38 need a bad-region mask like GRCh37's chr2/chr21?

Two independent tests on high_cov/ibd_all.ibd.gz:
  (a) genetic-map discontinuities -- big cM jumps over tiny bp windows, which
      turn any local IBS sharing into a spurious multi-cM "segment";
  (b) segment start/end pileups -- a real IBD boundary is a recombination
      event and should never repeat at one coordinate across many pairs.
"""
import gzip, glob, os
import numpy as np
from collections import Counter

B = os.path.dirname(os.path.abspath(__file__))

print("=== (a) genetic-map discontinuities (cM jump per 100 kb) ===")
worst = []
for c in list(range(1, 23)):
    d = np.loadtxt(f"{B}/genmap/plink.chr{c}.GRCh38.map", usecols=(2, 3))
    d = d[np.argsort(d[:, 1])]
    cm, bp = d[:, 0], d[:, 1]
    dcm, dbp = np.diff(cm), np.diff(bp)
    ok = dbp > 0
    rate = np.zeros(len(dcm)); rate[ok] = dcm[ok] / (dbp[ok] / 1e5)   # cM per 100kb
    i = int(np.argmax(rate))
    worst.append((rate[i], c, bp[i], dcm[i], dbp[i]))
worst.sort(reverse=True)
for r, c, pos, dc, db in worst[:6]:
    print(f"  chr{c:<3} {pos:>12,.0f}  {dc:6.2f} cM over {db:>9,.0f} bp  = {r:7.2f} cM/100kb")

print("\n=== (b) segment endpoint pileups ===")
starts, ends, tot = Counter(), Counter(), 0
with gzip.open(f"{B}/ibd_all.ibd.gz", "rt") as fh:
    for line in fh:
        c = line.split(); tot += 1
        starts[(c[4], c[5])] += 1; ends[(c[4], c[6])] += 1
print(f"  total segments: {tot:,}")
for nm, cnt in (("START", starts), ("END", ends)):
    top = cnt.most_common(5)
    print(f"  most common {nm} coordinates:")
    for (ch, pos), n in top:
        print(f"    {ch}:{int(pos):>12,}  {n:>7,} segments ({100*n/tot:.2f}%)")
