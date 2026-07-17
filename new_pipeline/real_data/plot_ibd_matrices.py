"""Plot the pairwise IBD spectra saved by generate_stan_data.py.

Produces TWO separate figures from one ``<prefix>.npz`` / ``<prefix>.json``, each
an n_pop x n_pop grid whose (i,j) cell is the length spectrum of population pair
(i,j) -- value vs. segment length (bin), on log-log axes:

  1. <prefix>_ibd_sharing.png : IBD *sharing* spectrum  (ibd_mean per length bin)
  2. <prefix>_ibd_count.png   : IBD *segment-count* spectrum (ibd_count per bin)

Usage
-----
    python plot_ibd_matrices.py <prefix> [--density] [--out-dir DIR]

``--density`` divides the sharing spectrum by bin width (IBD fraction / cM).
"""
import os
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load(prefix):
    d = np.load(prefix + ".npz")
    with open(prefix + ".json") as f:
        meta = json.load(f)
    return d, meta


def _spectrum_grid(stack, pops, bins, title, ylabel, out, logy=True,
                   width=None):
    """n_pop x n_pop grid; cell (i,j) is stack[:, i, j] vs bin midpoint."""
    P = len(pops)
    lo = bins[:, 0]
    hi = bins[:, 1]
    mid = np.sqrt(lo * hi)
    fig, axes = plt.subplots(P, P, figsize=(3.4 * P, 3.0 * P),
                             sharex=True, squeeze=False)
    for i in range(P):
        for j in range(P):
            ax = axes[i][j]
            if j < i:                       # symmetric -> hide lower triangle
                ax.axis("off")
                continue
            y = stack[:, i, j].astype(float)
            if width is not None:
                y = y / width
            nz = y > 0
            ax.plot(mid, y, "-", color="#c0392b", lw=1.4, alpha=0.6, zorder=1)
            ax.plot(mid[nz], y[nz], "o", color="k", ms=4, zorder=2)
            ax.set_xscale("log")
            if logy:
                ax.set_yscale("log")
            ax.set_title(f"{pops[i]}-{pops[j]}", fontsize=11)
            ax.grid(True, which="both", alpha=0.15)
            if i == P - 1 or j == i:
                ax.set_xlabel("segment length (cM)")
            if j == i:
                ax.set_ylabel(ylabel)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="Output prefix (without .npz/.json)")
    ap.add_argument("--density", action="store_true",
                    help="Divide the sharing spectrum by bin width (fraction/cM).")
    ap.add_argument("--out-dir", default=None,
                    help="Directory for the PNGs (default: alongside the prefix).")
    args = ap.parse_args()

    d, meta = _load(args.prefix)
    pops = meta["pop_order"]
    bins = np.array(meta["bins"], dtype=float)
    width = (bins[:, 1] - bins[:, 0]) if args.density else None
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.prefix))
    base = os.path.basename(args.prefix)

    if "ibd_mean" not in d.files:
        raise KeyError(f"{args.prefix}.npz has no 'ibd_mean' array.")
    ylab = "IBD fraction / cM" if args.density else "IBD fraction"
    _spectrum_grid(d["ibd_mean"], pops, bins,
                   f"Pairwise IBD sharing spectrum -- {base}",
                   ylab, os.path.join(out_dir, f"{base}_ibd_sharing.png"),
                   width=width)

    if "ibd_count" not in d.files:
        raise KeyError(f"{args.prefix}.npz has no 'ibd_count' array "
                       f"(regenerate with the current generate_stan_data.py).")
    _spectrum_grid(d["ibd_count"], pops, bins,
                   f"Pairwise IBD segment-count spectrum -- {base}",
                   "number of IBD segments",
                   os.path.join(out_dir, f"{base}_ibd_count.png"))


if __name__ == "__main__":
    main()
