"""Mask = HapNe chromosome arms (positive) + IBDNe's extreme-locus rule (negative).

(A) HapNe arms, GRCh38 (Fournier et al. 2023, Nat Commun 14:7107,
    src/hapne/files/regions_grch38.txt): 39 arms.  Everything outside them --
    centromeres, acrocentric short arms, telomere caps -- is excluded by
    construction rather than by any data-derived criterion.

(B) IBDNe's rule (Browning & Browning 2015, AJHG 97:404-418; documented at
    faculty.washington.edu/browning/ibdne.html):

        "Any genomic region that is less than 50 centiMorgans in length and
         bounded on each side by a chromosome boundary or a locus with an
         extreme number of IBD segments is excluded."

    IBDNe does not publish the constant behind "extreme", so it is chosen here
    from the data and the choice is reported with a sensitivity sweep.

Everything is done in GENETIC coordinates.  Under a neutral coalescent the number
of segments covering a locus is flat per cM, not per Mb, so a per-Mb density is
confounded by local recombination rate -- that confound is what made our earlier
criterion (d) hard to interpret.

Counting uses WITHIN-SUBPOPULATION pairs only, which is the setting IBDNe is
defined for (one population at a time); cross-population pairs carry almost no
true IBD and would change what "extreme" means.
"""
import os, gzip, argparse
import numpy as np

B = os.path.dirname(os.path.abspath(__file__))
HC = os.path.dirname(B)

ap = argparse.ArgumentParser()
ap.add_argument("--arms", default=f"{B}/hapne_arms_grch38.txt")
ap.add_argument("--ibd", default=f"{HC}/ibd_all.ibd.gz")
ap.add_argument("--labels", default=f"{HC}/samples/labels_subpop.txt")
ap.add_argument("--genmap-dir", default=f"{HC}/genmap")
ap.add_argument("--min-cm", type=float, default=2.0)
ap.add_argument("--grid", type=float, default=0.1, help="locus grid, cM")
ap.add_argument("--extreme", type=float, default=None,
                help="flag a locus above this multiple of the arm-wise median "
                     "coverage; if omitted, only the sweep is printed")
ap.add_argument("--fill-cm", type=float, default=50.0,
                help="IBDNe's rule: drop any stretch shorter than this bounded "
                     "on both sides by a flagged locus or an arm boundary")
ap.add_argument("--out", default=f"{B}/bad_regions_arm_ibdne_grch38.bed")
a = ap.parse_args()

# ---------------- genetic map ----------------
GM = {}
def gmap(ch):
    c = ch[3:] if ch.startswith("chr") else ch
    if c not in GM:
        d = np.loadtxt(f"{a.genmap_dir}/plink.chr{c}.GRCh38.map", usecols=(2, 3))
        GM[c] = d[np.argsort(d[:, 1])]
    return GM[c]

def bp2cm(ch, p):
    d = gmap(ch)
    return np.interp(p, d[:, 1], d[:, 0])

def cm2bp(ch, x):
    d = gmap(ch)
    return np.interp(x, d[:, 0], d[:, 1])

# ---------------- (A) arms ----------------
arms = []
for line in open(a.arms):
    if line.startswith("#") or line.startswith("CHR"):
        continue
    c = line.split()
    ch = f"chr{c[0]}"
    arms.append((ch, int(c[1]), int(c[2])))
arms.sort(key=lambda r: (int(r[0][3:]), r[1]))
arm_cm = {(ch, s, e): (bp2cm(ch, s), bp2cm(ch, e)) for ch, s, e in arms}
print(f"[A] {len(arms)} HapNe arms, "
      f"{sum(e - s for (_, _, _), (s, e) in arm_cm.items()):.0f} cM inside arms")

# ---------------- coverage on a cM grid, within-population pairs ----------------
sub = {}
for line in open(a.labels):
    if line.strip():
        s, p = line.split()[:2]
        sub[s] = p

grids = {}
for k, (s_cm, e_cm) in arm_cm.items():
    n = max(2, int(np.ceil((e_cm - s_cm) / a.grid)) + 1)
    grids[k] = np.zeros(n + 1)          # +1/-1 difference array

n_used = n_out = 0
with gzip.open(a.ibd, "rt") as fh:
    for line in fh:
        f = line.split()
        p1, p2 = sub.get(f[0]), sub.get(f[2])
        if p1 is None or p1 != p2 or float(f[7]) < a.min_cm:
            continue
        ch, s, e = f[4], int(f[5]), int(f[6])
        placed = False
        for key in arm_cm:
            if key[0] != ch or e <= key[1] or s >= key[2]:
                continue
            s_cm, e_cm = arm_cm[key]
            x0 = max(bp2cm(ch, max(s, key[1])), s_cm)
            x1 = min(bp2cm(ch, min(e, key[2])), e_cm)
            i0 = int((x0 - s_cm) / a.grid)
            i1 = int((x1 - s_cm) / a.grid)
            g = grids[key]
            i0 = max(0, min(i0, len(g) - 2)); i1 = max(i0 + 1, min(i1, len(g) - 1))
            g[i0] += 1; g[i1] -= 1
            placed = True
        n_used += placed
        n_out += (not placed)
print(f"[cov] {n_used:,} within-population segments >= {a.min_cm} cM placed inside arms "
      f"({n_out:,} fell entirely outside)")

cov = {k: np.cumsum(v)[:-1] for k, v in grids.items()}
allcov = np.concatenate(list(cov.values()))
med = float(np.median(allcov[allcov > 0])) if (allcov > 0).any() else 0.0
print(f"[cov] genome-wide median coverage over occupied loci = {med:.0f} segments")
qs = [50, 75, 90, 95, 99, 99.5, 99.9, 100]
print("      coverage quantiles (x median): "
      + "  ".join(f"p{q}={np.percentile(allcov, q)/med:.1f}" for q in qs))

# ---------------- threshold sweep ----------------
def flag_and_fill(thresh):
    """Flag loci above thresh*median, then apply IBDNe's <fill-cm bounded rule."""
    out = []
    for key in cov:
        ch, bs, be = key
        s_cm, e_cm = arm_cm[key]
        c = cov[key]
        bad = c > thresh * med
        # boundaries of the run: arm start, flagged runs, arm end
        idx = np.flatnonzero(bad)
        # blocks of contiguous flagged loci
        blocks = []
        if len(idx):
            st = idx[0]; pv = idx[0]
            for i in idx[1:]:
                if i == pv + 1:
                    pv = i
                else:
                    blocks.append((st, pv)); st = i; pv = i
            blocks.append((st, pv))
        # anchors = arm boundaries + flagged blocks
        anchors = [(-1, -1)] + blocks + [(len(c), len(c))]
        keep = []
        for (a0, a1), (b0, b1) in zip(anchors[:-1], anchors[1:]):
            gap_lo, gap_hi = a1 + 1, b0            # unflagged stretch between them
            if gap_hi <= gap_lo:
                continue
            keep.append((gap_lo, gap_hi, (gap_hi - gap_lo) * a.grid))
        # everything not kept, plus kept stretches shorter than fill-cm
        drop = [(b0, b1 + 1) for b0, b1 in blocks]
        drop += [(lo, hi) for lo, hi, L in keep if L < a.fill_cm]
        for lo, hi in drop:
            x0 = s_cm + lo * a.grid
            x1 = s_cm + min(hi, len(c)) * a.grid
            out.append((ch, int(cm2bp(ch, x0)), int(cm2bp(ch, x1)), x1 - x0))
    return out

print(f"\n[sweep] IBDNe rule: flag loci > T x median, then drop unflagged stretches "
      f"< {a.fill_cm:g} cM between anchors")
print(f"  {'T':>6} {'loci flagged':>13} {'regions':>8} {'cM dropped':>11} {'% of arms':>10}")
arm_total = sum(e - s for (s, e) in arm_cm.values())
for T in (2, 3, 4, 5, 6, 8, 10, 15, 20):
    regs = flag_and_fill(T)
    cm_drop = sum(r[3] for r in regs)
    nflag = int(sum((c > T * med).sum() for c in cov.values()))
    print(f"  {T:>6.0f} {nflag:>13,} {len(regs):>8,} {cm_drop:>11.1f} "
          f"{100*cm_drop/arm_total:>9.1f}%")

if a.extreme is None:
    print("\nNo --extreme given; sweep only.  Re-run with a threshold to write the BED.")
    raise SystemExit

# ---------------- final mask ----------------
regs = flag_and_fill(a.extreme)
# outside-arm exclusions
CHRLEN = {f"chr{c}": int(gmap(f"chr{c}")[-1, 1]) + 5_000_000 for c in range(1, 23)}
outside = []
by_chr = {}
for ch, s, e in arms:
    by_chr.setdefault(ch, []).append((s, e))
for ch in sorted(by_chr, key=lambda x: int(x[3:])):
    prev = 0
    for s, e in sorted(by_chr[ch]):
        if s > prev:
            outside.append((ch, prev, s))
        prev = e
    if CHRLEN[ch] > prev:
        outside.append((ch, prev, CHRLEN[ch]))

allr = [(ch, s, e, "outside-arm") for ch, s, e in outside] \
     + [(ch, s, e, "extreme-ibd") for ch, s, e, _ in regs]
merged = []
for ch, s, e, tag in sorted(allr, key=lambda r: (int(r[0][3:]), r[1])):
    if merged and merged[-1][0] == ch and s <= merged[-1][2]:
        merged[-1][2] = max(merged[-1][2], e)
        if tag not in merged[-1][3]:
            merged[-1][3] += "+" + tag
    else:
        merged.append([ch, s, e, tag])

with open(a.out, "w") as f:
    f.write("# IBD mask, GRCh38 = HapNe arms + IBDNe extreme-locus rule\n")
    f.write("# outside-arm : complement of the 39 HapNe chromosome arms "
            "(Fournier et al. 2023, regions_grch38.txt)\n")
    f.write(f"# extreme-ibd : IBDNe rule (Browning & Browning 2015) -- locus coverage "
            f"> {a.extreme}x the median ({med:.0f} segments), then any unflagged stretch "
            f"< {a.fill_cm:g} cM between anchors also dropped\n")
    for ch, s, e, tag in merged:
        f.write(f"{ch}\t{s}\t{e}\t{tag}\n")
tot = sum(bp2cm(c, e) - bp2cm(c, s) for c, s, e, _ in merged)
print(f"\n[save] {a.out}: {len(merged)} regions, {tot:.1f} cM "
      f"({100*tot/3372:.1f}% of the genetic map)")
for ch, s, e, tag in merged:
    if tag != "outside-arm":
        print(f"   {ch}:{s/1e6:.1f}-{e/1e6:.1f} Mb  "
              f"{bp2cm(ch,e)-bp2cm(ch,s):6.2f} cM  {tag}")
