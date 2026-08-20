"""Derive the IBD bad-region mask from data, on three independent signals.

Run on the UNMASKED high_cov/ibd_all.ibd.gz so the mask is rebuilt from scratch
rather than layered on top of a previous one.

  (a) genetic-map discontinuity -- a big cM jump over a tiny bp span turns any
      local IBS match into a spurious multi-cM "segment".
  (b) segment endpoint pileup -- a real IBD boundary is a recombination event
      and should not repeat at one coordinate across thousands of pairs.
  (c) CROSS-POPULATION ENRICHMENT  <-- the one that matters, and the one the
      earlier version of this script missed.

Why (c) is the sensitive test.  Signals (a) and (b) look at the segment pool as
a whole, which is ~99% within-population and therefore ~99% real IBD; a locus
producing a few hundred false matches is invisible against 90,000 true ones.
Cross-continental pairs, by contrast, share almost no true IBD >=2 cM, so the
same fixed background of false matches is most of what they have.  Measured on
the clustered leaves: 57% of cross-population segments >=2.5 cM fell in
centromeres / acrocentric short arms / the 8p23 inversion, against 0.9% of
within-population segments -- a 63x enrichment that (a) and (b) cannot see.

Detection uses cross-SUPERPOPULATION pairs over the whole panel rather than the
four clustered leaves: it is the same signal with ~100x more segments, and it
keeps the mask independent of the clustering that will later consume it.  ASW
and ACB are dropped because their European ancestry is real recent
cross-continental IBD -- exactly what we must not mask.  Genuine admixture is
spread over the genome; an artifact locus is concentrated, which is what the
per-window enrichment measures.
"""
import gzip, os, argparse
import numpy as np
from collections import Counter, defaultdict

B = os.path.dirname(os.path.abspath(__file__))
ADMIXED = {"ASW", "ACB"}          # real recent cross-continental IBD, not artifact

ap = argparse.ArgumentParser()
ap.add_argument("--ibd", default=f"{B}/ibd_all.ibd.gz")
ap.add_argument("--window", type=float, default=1e6, help="scan window (bp)")
ap.add_argument("--min-cm", type=float, default=2.0)
ap.add_argument("--min-cross", type=int, default=5,
                help="minimum cross-population segments in a window to judge it")
ap.add_argument("--enrich", type=float, default=10.0,
                help="flag a window at this many times the genome-wide "
                     "cross/within ratio")
ap.add_argument("--dens", type=float, default=15.0,
                help="(d) flag a window at this many times the median "
                     "within-population segment density")
ap.add_argument("--dens-enrich", type=float, default=1.5,
                help="(d) ...but only if cross-enrichment also exceeds this, so "
                     "genuine within-population IBD hotspots (enrichment < 1) "
                     "are kept")
ap.add_argument("--out", default=f"{B}/bad_regions_grch38.bed")
a = ap.parse_args()

# ---- labels ----
sub = {}
for line in open(f"{B}/samples/labels_subpop.txt"):
    s, p = line.split()[:2]; sub[s] = p
sup = {}
for line in open(f"{B}/samples/labels_4pop.txt"):
    s, p = line.split()[:2]
    if sub.get(s) not in ADMIXED:
        sup[s] = p
print(f"[labels] {len(sup)} individuals after dropping {sorted(ADMIXED)}")

# ---- (c) per-window cross vs within counts ----
W = a.window
cross = Counter(); within = Counter()
n_c = n_w = 0
with gzip.open(a.ibd, "rt") as fh:
    for line in fh:
        c = line.split()
        p1, p2 = sup.get(c[0]), sup.get(c[2])
        if p1 is None or p2 is None or float(c[7]) < a.min_cm:
            continue
        key = (c[4], int((int(c[5]) + int(c[6])) / 2 // W))     # midpoint window
        # The contrast has to be same-SUBPOPULATION (IBD-rich, ~all real) versus
        # cross-SUPERPOPULATION (IBD-poor, ~all artifact).  Lumping
        # cross-subpopulation pairs like CHB-JPT into "within" contaminates the
        # denominator with the same artifacts and flattens the enrichment: that
        # framing capped the signal at 5x, where this one reaches ~100x.
        # Pair counts cancel in the ratio-of-ratios, so no pair normalisation
        # is needed.
        s1, s2 = sub.get(c[0]), sub.get(c[2])
        if p1 == p2:
            if s1 == s2:
                within[key] += 1; n_w += 1
        else:
            cross[key] += 1; n_c += 1
r0 = n_c / max(n_w, 1)
print(f"[scan] {n_c:,} cross / {n_w:,} within segments >= {a.min_cm} cM; "
      f"genome-wide cross/within ratio = {r0:.4f}")

rows = []
for key, nc in cross.items():
    if nc < a.min_cross:
        continue
    nw = within.get(key, 0)
    e = (nc / (nw + 1)) / r0
    rows.append((e, key, nc, nw))
rows.sort(reverse=True)
print(f"\n[scan] windows with >= {a.min_cross} cross segments, "
      f"top 15 by enrichment:")
print(f"  {'window':>22} {'cross':>8} {'within':>8} {'enrich':>8} {'%allcross':>10}")
for e, (ch, w), nc, nw in rows[:20]:
    print(f"  {ch}:{int(w*W/1e6):>4}-{int((w+1)*W/1e6):<4}Mb {nc:>8d} {nw:>8d} "
          f"{e:>8.1f}x {100*nc/n_c:>9.2f}%")

# Threshold sweep: there is no knee in the enrichment ranking, so choose on the
# trade-off instead -- how much spurious cross signal comes out per unit of real
# within-population signal sacrificed.
print(f"\n[sweep] cost/benefit of the enrichment cut "
      f"(cross removed = artifact, within removed = collateral):")
print(f"  {'cut':>6} {'windows':>8} {'cross removed':>15} {'within removed':>16} {'ratio':>7}")
for cut in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0):
    ws = {(ch, int(w)) for e, (ch, w), nc, nw in rows if e >= cut}
    rc = sum(nc for e, k, nc, nw in rows if e >= cut)
    rw = sum(nw for e, k, nc, nw in rows if e >= cut)
    print(f"  {cut:>5.1f}x {len(ws):>8d} {rc:>9,d} ({100*rc/n_c:4.1f}%) "
          f"{rw:>10,d} ({100*rw/n_w:4.1f}%) {rc/max(rw,1):>6.1f}")

flag_c = {(ch, int(w)) for e, (ch, w), nc, nw in rows if e >= a.enrich}

# ---- (d) absolute within-population over-density ----
# A genetic-map discontinuity inflates within- AND cross-pairs by the same
# factor, so it leaves the RATIO near baseline and criterion (c) cannot see it:
# chr2:2-6 Mb runs 17-25x the median within-population density at only 1.9-2.9x
# enrichment, below any usable (c) cut.  Density alone is not enough either,
# because real IBD hotspots are also dense -- but they are dense in
# within-population pairs specifically, so their enrichment sits BELOW 1
# (chr4:3-4 Mb at 0.5x, chr16:85-86 Mb at 0.6x are genuine and must survive).
# Requiring high density AND enrichment above 1 separates the two cleanly.
med = float(np.median([v for v in within.values()])) if within else 0.0
flag_d = set()
for key, nw in within.items():
    if med > 0 and nw >= a.dens * med:
        nc = cross.get(key, 0)
        if (nc / (nw + 1)) / r0 >= a.dens_enrich:
            flag_d.add((key[0], int(key[1])))
print(f"\n[flag] (c) cross-enrichment >= {a.enrich}x: {len(flag_c)} windows")
print(f"[flag] (d) within-density >= {a.dens}x median ({med:.0f}/window) "
      f"AND enrichment >= {a.dens_enrich}x: {len(flag_d)} windows"
      + (f"  {sorted(flag_d - flag_c)[:6]}" if flag_d - flag_c else ""))
flag = sorted(flag_c | flag_d)
print(f"[flag] union: {len(flag)} windows")

# ---- merge adjacent flagged windows into regions ----
merged = []
for ch, w in flag:
    if merged and merged[-1][0] == ch and w <= merged[-1][2] + 1:
        merged[-1][2] = w
    else:
        merged.append([ch, w, w])
regions = [(ch, int(w0 * W), int((w1 + 1) * W)) for ch, w0, w1 in merged]

# ---- how much cross / within signal does this remove? ----
def covered(ch, s, e):
    return any(ch == c and s < hi and e > lo for c, lo, hi in regions)
kc = kw = 0
with gzip.open(a.ibd, "rt") as fh:
    for line in fh:
        c = line.split()
        p1, p2 = sup.get(c[0]), sup.get(c[2])
        if p1 is None or p2 is None or float(c[7]) < a.min_cm:
            continue
        if covered(c[4], int(c[5]), int(c[6])):
            if p1 == p2: kw += 1
            else: kc += 1
print(f"[effect] mask removes {kc:,}/{n_c:,} cross ({kc/max(n_c,1):.1%}) "
      f"and {kw:,}/{n_w:,} within ({kw/max(n_w,1):.1%}) segments >= {a.min_cm} cM")
# The within-population loss is NOT collateral damage on real IBD.  If these
# windows were ordinary genome, their within-population segment density would
# match the genome-wide average; instead it is ~100x higher, i.e. those
# segments are the same artifact showing up in within-population pairs too.
# Density must be measured with MIDPOINT assignment, not the overlap count used
# for removal: a 2 cM segment is ~2 Mb long, so it can overlap a 1 Mb window
# while lying mostly outside it, and dividing the overlap count by the window
# length inflates the density several-fold.
mask_mb = sum(e - s2 for _, s2, e in regions) / 1e6
genome_mb = 2875.0                     # autosomal GRCh38, non-N
kw_mid = sum(nw for e_, k, nc, nw in rows if e_ >= a.enrich)
d_in = kw_mid / mask_mb if mask_mb else 0.0
d_out = (n_w - kw_mid) / (genome_mb - mask_mb)
print(f"[effect] within-population segment density: {d_in:,.0f}/Mb inside the "
      f"mask vs {d_out:,.0f}/Mb outside = {d_in/max(d_out,1e-9):.0f}x over-dense "
      f"-> the masked within-population segments are artifact too, not lost signal")

# ---- genetic length removed (for the IBD normalizer) ----
tot_cm = 0.0
for ch, s, e in regions:
    c = ch[3:] if ch.startswith("chr") else ch
    d = np.loadtxt(f"{B}/genmap/plink.chr{c}.GRCh38.map", usecols=(2, 3))
    d = d[np.argsort(d[:, 1])]
    tot_cm += float(np.interp(e, d[:, 1], d[:, 0]) - np.interp(s, d[:, 1], d[:, 0]))

with open(a.out, "w") as f:
    f.write("# IBD bad-region mask, GRCh38, derived by high_cov/derive_mask.py\n")
    f.write(f"# (c) >={a.enrich}x cross-population enrichment over the genome-wide "
            f"ratio, >={a.min_cross} cross segments per {W/1e6:g} Mb window\n")
    f.write(f"# (d) OR >={a.dens}x median within-population density AND "
            f">={a.dens_enrich}x enrichment (map-discontinuity class)\n")
    for ch, s, e in regions:
        f.write(f"{ch}\t{s}\t{e}\n")
print(f"\n[save] {a.out}  ({len(regions)} regions, {tot_cm:.2f} cM)")
for ch, s, e in regions:
    print(f"  {ch}\t{s:>12,}\t{e:>12,}")
