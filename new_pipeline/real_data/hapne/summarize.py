"""Collect every HapNe-IBD fit under hapne/ into one table and one figure.

Reported alongside Ne: the two things that predict whether a HapNe fit is
trustworthy at all.

  seg/pair  IBD segments per within-population pair.  Established earlier on 20
            fits: r = -0.98 against log10 Ne(100), with a clean split at 2.5 --
            above it Ne(100) lands in the 13.5k-28.6k textbook band, below it in
            70.7k-585k.  A homogeneous ~100-sample draw from a large outbred
            population simply cannot constrain recent Ne, and no amount of
            modelling rescues that.

  arms      How many of the 39 chromosome arms carried any usable segment.
            HapNe estimates its per-bin variance from the spread ACROSS arms
            (supp. Eq 25-26), so a fit resting on few arms has an unreliable
            covariance even where the point estimate looks fine.

Do NOT judge these curves by shape.  Fold-change and monotonicity both rank the
RAILED fits highest -- a curve that rails to 1.8e7 and then falls smoothly scores
83x fold-change and 1.00 monotonicity.  Judge by the LEVEL: Ne at generation 100
against the ~1e4 ancestral human value.
"""
import os, gzip, glob, json
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = os.path.dirname(os.path.abspath(__file__))
HC = os.path.join(os.path.dirname(B), "high_cov")

sub, sup = {}, {}
for line in open(f"{HC}/samples/labels_subpop.txt"):
    s, p = line.split()[:2]; sub[s] = p
for line in open(f"{HC}/samples/labels_4pop.txt"):
    s, p = line.split()[:2]; sup[s] = p
size = Counter(sub.values()); size.update()
size_sup = Counter(sup.values())

seg = Counter()
with gzip.open(f"{HC}/ibd_all_masked.ibd.gz", "rt") as fh:
    for line in fh:
        f = line.split()
        a, b = f[0], f[2]
        if sub.get(a) and sub.get(a) == sub.get(b):
            seg[sub[a]] += 1
        if sup.get(a) and sup.get(a) == sup.get(b):
            seg[sup[a]] += 1

rows = []
for d in sorted(glob.glob(f"{B}/*/HapNe/hapne.csv")):
    name = os.path.basename(os.path.dirname(os.path.dirname(d)))
    r = pd.read_csv(d)
    n = size.get(name) or size_sup.get(name)
    if not n:
        continue
    npair = n * (n - 1) / 2
    arms = len(glob.glob(f"{B}/{name}/IBD/*.ibd.gz"))
    get = lambda g: (float(r["Q0.5"][r.TIME == g].iloc[0]) if (r.TIME == g).any()
                     else float("nan"))
    rows.append({"pop": name, "kind": "super" if name in size_sup and name not in size else "sub",
                 "n": n, "segs": seg[name], "seg_per_pair": seg[name] / npair,
                 "arms": arms, "Ne10": get(10), "Ne25": get(25),
                 "Ne50": get(50), "Ne100": get(100)})

df = pd.DataFrame(rows).sort_values(["kind", "seg_per_pair"], ascending=[True, False])
df.to_csv(f"{B}/hapne_summary.csv", index=False)

print(f"{'pop':>6} {'kind':>6} {'n':>5} {'segs':>8} {'seg/pair':>9} {'arms':>5} "
      f"{'Ne(10)':>10} {'Ne(25)':>10} {'Ne(50)':>10} {'Ne(100)':>10}  verdict")
for _, x in df.iterrows():
    ok = 8e3 <= x.Ne100 <= 4e4
    v = "sane" if ok else ("RAILED" if x.Ne100 > 4e4 else "low")
    print(f"{x['pop']:>6} {x['kind']:>6} {x.n:>5.0f} {x.segs:>8,.0f} "
          f"{x.seg_per_pair:>9.3f} {x.arms:>5.0f} "
          + "".join(f"{x[c]:>10,.0f}" for c in ("Ne10", "Ne25", "Ne50", "Ne100"))
          + f"  {v}")

d = df.dropna(subset=["Ne100"])
if len(d) > 3:
    r = np.corrcoef(np.log10(d.seg_per_pair.clip(1e-3)), np.log10(d.Ne100.clip(1)))[0, 1]
    print(f"\ncorr(log10 seg/pair, log10 Ne(100)) = {r:+.2f}   "
          f"[established value on the previous mask: -0.98]")
    hi, lo = d[d.seg_per_pair >= 2.5], d[d.seg_per_pair < 2.5]
    if len(hi) and len(lo):
        print(f"   seg/pair >= 2.5 (n={len(hi):2d}): Ne(100) "
              f"{hi.Ne100.min():>9,.0f} - {hi.Ne100.max():>9,.0f}")
        print(f"   seg/pair <  2.5 (n={len(lo):2d}): Ne(100) "
              f"{lo.Ne100.min():>9,.0f} - {lo.Ne100.max():>9,.0f}")

fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
for _, x in df.iterrows():
    f = f"{B}/{x['pop']}/HapNe/hapne.csv"
    r = pd.read_csv(f)
    k = 0 if x["kind"] == "sub" else 1
    ax[k].plot(r.TIME, r["Q0.5"], lw=1.4, label=f"{x['pop']} ({x.seg_per_pair:.1f})")
for k, t in enumerate(["1000G subpopulations", "superpopulations"]):
    ax[k].set_yscale("log"); ax[k].set_xlabel("generations ago")
    ax[k].set_ylabel("Ne (diploid)"); ax[k].set_title(f"HapNe-IBD, {t}", fontsize=11)
    ax[k].axhspan(8e3, 4e4, color="grey", alpha=.15, zorder=0)
    ax[k].legend(fontsize=7, ncol=2)
fig.suptitle("HapNe-IBD after the referenced mask (arms + IBDNe extreme loci + "
             "cross-population enrichment); shaded = plausible ancestral Ne", fontsize=12)
fig.tight_layout()
fig.savefig(f"{B}/hapne_summary.png", dpi=130)
print(f"\n[save] {B}/hapne_summary.csv  and  hapne_summary.png")
