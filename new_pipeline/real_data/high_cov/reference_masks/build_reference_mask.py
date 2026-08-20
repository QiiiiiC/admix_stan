"""Build an IBD analysis mask from PUBLISHED region definitions, and audit our
data-derived mask against it.

Two independent published sources, which do different jobs:

  (A) HapNe chromosome ARMS (Fournier et al. 2023, Nat Commun 14:7107;
      src/hapne/files/regions_grch38.txt).  A POSITIVE definition of the
      analysable genome: 39 arms, centromeres and acrocentric short arms dropped
      by construction.  This is the convention used by the very method we
      benchmark our Ne against, so adopting it makes the comparison like-for-like
      instead of us inventing our own footprint.

  (B) Li et al. 2014 excess-IBD regions (PLoS Genet 10(1):e1004144, Table 3).
      A NEGATIVE list: 14 regions >= 5 cM where observed/expected IBD length
      exceeds 4x, measured with GERMLINE, ISCA and fastIBD across European, East
      Asian and Mexican American controls.  These are inside the arms and would
      otherwise survive (A).  Published in hg19, so lifted here.

Note on what Browning & Browning (2015, AJHG 97:404-418) actually provide: NOT a
static region list.  IBDNe excludes dynamically -- "any genomic region that is
less than 50 centiMorgans in length and bounded on each side by a chromosome
boundary or a locus with an extreme number of IBD segments is excluded" -- plus
the mincm floor (default 2 cM) and the recommendation to run merge-ibd-segments
first.  So it cannot be used as a reference BED; its analogue here is (B).

Output: bad_regions_reference_grch38.bed, the complement of (A) intersected with
the autosomes, unioned with lifted (B).
"""
import os, gzip, argparse
import numpy as np

B = os.path.dirname(os.path.abspath(__file__))
HC = os.path.dirname(B)

ap = argparse.ArgumentParser()
ap.add_argument("--arms", default=f"{B}/hapne_arms_grch38.txt")
ap.add_argument("--li", default=f"{B}/li2014_excess_ibd_hg19.bed")
ap.add_argument("--chain", default=f"{B}/hg19ToHg38.over.chain.gz")
ap.add_argument("--ours", default=f"{HC}/bad_regions_grch38.bed")
ap.add_argument("--genmap-dir", default=f"{HC}/genmap")
ap.add_argument("--out", default=f"{B}/bad_regions_reference_grch38.bed")
a = ap.parse_args()

# ---------- (A) arms -> the analysable footprint, per chromosome ----------
arms = {}
for line in open(a.arms):
    if line.startswith("#") or line.startswith("CHR"):
        continue
    c = line.split()
    arms.setdefault(f"chr{c[0]}", []).append((int(c[1]), int(c[2])))
for k in arms:
    arms[k].sort()
print(f"[A] HapNe arms: {sum(len(v) for v in arms.values())} regions over "
      f"{len(arms)} autosomes")

# ---------- (B) Li et al., lifted hg19 -> GRCh38 ----------
from pyliftover import LiftOver
lo = LiftOver(a.chain)


def lift(ch, pos):
    r = lo.convert_coordinate(ch, pos)
    return r[0][1] if r else None


li38, li_fail = [], []
for line in open(a.li):
    if line.startswith("#") or not line.strip():
        continue
    c = line.split()
    ch, s, e = c[0], int(c[1]), int(c[2])
    s2, e2 = lift(ch, s), lift(ch, e)
    if s2 is None or e2 is None or e2 <= s2:
        # Endpoints inside sequence that changed between builds do not lift.
        # Walk inward in 100 kb steps rather than silently dropping the region.
        step = 100_000
        for off in range(0, 3_000_000, step):
            if s2 is None:
                s2 = lift(ch, s + off)
            if e2 is None:
                e2 = lift(ch, e - off)
            if s2 is not None and e2 is not None and e2 > s2:
                break
    if s2 is None or e2 is None or e2 <= s2:
        li_fail.append((ch, s, e))
        continue
    li38.append((ch, s2, e2))
print(f"[B] Li et al. 2014: {len(li38)}/14 regions lifted to GRCh38"
      + (f"; FAILED {li_fail}" if li_fail else ""))


def genmap(ch):
    c = ch[3:] if ch.startswith("chr") else ch
    d = np.loadtxt(f"{a.genmap_dir}/plink.chr{c}.GRCh38.map", usecols=(2, 3))
    return d[np.argsort(d[:, 1])]


GM = {}
def cm(ch, s, e):
    if ch not in GM:
        GM[ch] = genmap(ch)
    d = GM[ch]
    return float(np.interp(e, d[:, 1], d[:, 0]) - np.interp(s, d[:, 1], d[:, 0]))


# ---------- complement of the arms = everything (A) already excludes ----------
CHRLEN = {f"chr{c}": int(genmap(f"chr{c}")[-1, 1]) + 5_000_000 for c in range(1, 23)}
excl = []
for ch in sorted(arms, key=lambda x: int(x[3:])):
    prev = 0
    for s, e in arms[ch]:
        if s > prev:
            excl.append((ch, prev, s, "outside-arm"))
        prev = e
    if CHRLEN[ch] > prev:
        excl.append((ch, prev, CHRLEN[ch], "outside-arm"))
for ch, s, e in li38:
    excl.append((ch, s, e, "li2014"))


def merge(regs):
    out = []
    for ch, s, e, tag in sorted(regs, key=lambda r: (int(r[0][3:]), r[1])):
        if out and out[-1][0] == ch and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
            if tag not in out[-1][3]:
                out[-1][3] += "+" + tag
        else:
            out.append([ch, s, e, tag])
    return out


merged = merge(excl)
tot_cm = sum(cm(c, s, e) for c, s, e, _ in merged)
with open(a.out, "w") as f:
    f.write("# IBD bad-region mask from PUBLISHED sources, GRCh38.\n")
    f.write("# outside-arm : complement of HapNe's 39 chromosome arms "
            "(Fournier et al. 2023, regions_grch38.txt)\n")
    f.write("# li2014      : excess-IBD regions, Li et al. 2014 PLoS Genet "
            "10(1):e1004144 Table 3, lifted hg19->GRCh38\n")
    for ch, s, e, tag in merged:
        f.write(f"{ch}\t{s}\t{e}\t{tag}\n")
print(f"\n[save] {a.out}: {len(merged)} regions, {tot_cm:.1f} cM")

# ---------- audit our data-derived mask against the published one ----------
ours = []
for line in open(a.ours):
    if line.startswith("#") or not line.strip():
        continue
    c = line.split()
    ours.append((c[0], int(c[1]), int(c[2])))

def overlap(ch, s, e, regs):
    hits = [(rs, re_, tg) for rc, rs, re_, tg in regs if rc == ch and s < re_ and e > rs]
    if not hits:
        return 0.0, ""
    cov = sum(min(e, re_) - max(s, rs) for rs, re_, tg in hits)
    return cov / (e - s), "+".join(sorted({t for _, _, t in hits}))

print(f"\n{'our region':>22} {'cM':>6} {'covered by published':>21}  source")
n_cov = n_new = 0
for ch, s, e in sorted(ours, key=lambda r: (int(r[0][3:]), r[1])):
    frac, tag = overlap(ch, s, e, merged)
    lab = "NOVEL" if frac < 0.5 else f"{frac:.0%}"
    if frac >= 0.5:
        n_cov += 1
    else:
        n_new += 1
    print(f"{ch+':'+str(s//10**6)+'-'+str(e//10**6)+'Mb':>22} {cm(ch,s,e):>6.2f} "
          f"{lab:>21}  {tag or '-'}")
print(f"\n{n_cov}/{len(ours)} of our regions are corroborated by published sources; "
      f"{n_new} are novel.")
