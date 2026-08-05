"""Pooled superpopulation vs cluster-selected Ne, one panel per superpopulation.

Everything except the sample set is held fixed: same ibd_all_masked.ibd.gz, same
mask, same 39 arms, same HapNe settings.  So any difference between the curves
in a panel is attributable to who is in the population, not to the data or the
method.  Five series per panel: the pooled superpopulation, and the cluster
selected by each of the four clustering routes.

Two numbers per curve:

  seg/pair   IBD segments per within-population PAIR.  This, not sample size, is
             what a HapNe fit actually spends: two populations of ~100 people
             differ 30-fold in how much IBD they carry, and one that carries
             almost none cannot constrain the recent generations however
             homogeneous it is.
  Ne(100)    the deep end of the curve.  Used as the quality read-out because
             shape statistics do not work here -- fold change and monotonicity
             both SCORE THE RAILED FITS HIGHEST (EAS rails to 1.8e7 and then
             falls smoothly, giving 83x fold and monotonicity 1.00 while being
             obviously wrong).  The level is what separates them: every fit that
             lands near the ~1e4 ancestral human Ne is sane, every fit stuck at
             1e5-1e6 is not.
"""
import os
import glob
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POPS = ["AFR", "EAS", "EUR", "SAS"]
HN = os.path.join(B, "clustering", "hapne", "{p}_%s")
SERIES = [("pooled superpop",  os.path.join(B, "hapne", "{p}"), "#7f7f7f"),
          ("SNP kmeans",       HN % "snp",      "#1f77b4"),
          ("IBD louvain",      HN % "louvain",  "#c0392b"),
          ("IBD spectral",     HN % "ibd",      "#e08214"),
          ("IBD dbscan",       HN % "dbscan",   "#7b3294")]


def load(d):
    f = os.path.join(d, "HapNe", "hapne.csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def nsamples(d):
    cfg = os.path.join(d, "config.ini")
    if os.path.exists(cfg):
        for line in open(cfg):
            if line.startswith("nb_samples"):
                return int(line.split("=")[1])
    return None


def density(d, n):
    """(total segments, segments per within-population pair).

    Counts NON-BLANK lines.  Runs written before the run_pop.py newline fix have
    a blank line after every record, so a naive line count doubles them -- and
    since only some directories were rewritten, the error is not even a constant
    factor across the table.  HapNe itself is unaffected (its awk bins a blank
    line at 0, which it never prints), so only this diagnostic needed the guard.
    """
    tot = 0
    for f in glob.glob(os.path.join(d, "IBD", "*.ibd.gz")):
        with gzip.open(f, "rt") as fh:
            tot += sum(1 for ln in fh if ln.strip())
    npair = n * (n - 1) / 2 if n and n > 1 else np.nan
    return tot, (tot / npair if npair and npair > 0 else np.nan)


def summarise(d):
    r = load(d)
    if r is None:
        return None
    n = nsamples(d)
    nseg, dens = density(d, n)
    return dict(r=r, n=n, nseg=nseg, dens=dens,
                g1=float(r["Q0.5"].iloc[0]), g100=float(r["Q0.5"].iloc[-1]))


fig, axes = plt.subplots(1, len(POPS), figsize=(4.6 * len(POPS), 4.6),
                         sharey=True, squeeze=False)
rows = []
for c, p in enumerate(POPS):
    ax = axes[0][c]
    for name, tmpl, col in SERIES:
        s = summarise(tmpl.format(p=p))
        if s is None:
            continue
        rows.append((p, name, s))
        ax.plot(s["r"].TIME, s["r"]["Q0.5"], color=col, lw=2.0,
                label=f"{name}  n={s['n']}, {s['dens']:.1f} seg/pair", zorder=3)
        ax.fill_between(s["r"].TIME, s["r"]["Q0.025"], s["r"]["Q0.975"],
                        color=col, alpha=0.13, lw=0, zorder=2)
    # ~1e4-2e4 is where every well-constrained fit here lands at generation 100
    ax.axhspan(1e4, 2e4, color="#00a050", alpha=0.08, lw=0, zorder=0)
    ax.set_yscale("log")
    ax.set_xlabel("generations before present")
    ax.set_title(p, fontsize=12)
    ax.grid(alpha=0.25, which="both", lw=0.5)
    ax.legend(fontsize=8, loc="upper right")
axes[0][0].set_ylabel("effective population size $N_e$")
fig.suptitle(
    "HapNe-IBD: pooled 1000G superpopulation vs the cohesive cluster found by "
    "each clustering route\nidentical IBD data, mask and settings -- only the "
    "sample set differs;  green band = the ~1-2e4 ancestral $N_e$",
    fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.89])
out = os.path.join(B, "clustering", "hapne_cluster_comparison.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"[save] {out}")

allrows = rows
print(f"\n{'pop':4s} {'series':16s} {'n':>5s} {'segs':>8s} {'seg/pair':>9s} "
      f"{'Ne(1)':>12s} {'Ne(100)':>10s}")
for p, name, s in sorted(allrows, key=lambda t: -t[2]["dens"]):
    print(f"{p:4s} {name:16s} {s['n']:>5} {s['nseg']:8,d} {s['dens']:9.2f} "
          f"{s['g1']:12,.0f} {s['g100']:10,.0f}")

d = np.array([s["dens"] for _, _, s in allrows])
g = np.array([s["g100"] for _, _, s in allrows])
print(f"\nlog10(seg/pair) vs log10(Ne at gen 100):  r = "
      f"{np.corrcoef(np.log10(d), np.log10(g))[0, 1]:.2f}  over {len(d)} fits")
hi, lo = d >= 2.5, d < 2.5
print(f"  seg/pair >= 2.5 ({hi.sum()} fits): Ne(100) = "
      + ", ".join(f"{v:,.0f}" for v in np.sort(g[hi])))
print(f"  seg/pair <  2.5 ({lo.sum()} fits): Ne(100) = "
      f"{g[lo].min():,.0f} - {g[lo].max():,.0f}")
