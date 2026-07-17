"""Find a homogeneous sub-cluster within each superpopulation, two ways.

The coalescent model assumes any two individuals *within* a leaf population can
coalesce (panmixia).  That is violated when a 1000G superpopulation pools several
subpopulations (e.g. AFR = YRI/LWK/GWD/MSL/ESN + admixed ASW/ACB).  This script
detects that substructure and picks the dominant homogeneous cluster per
population, by two independent signals:

  (a) IBD:  cluster individuals on their pairwise IBD co-ancestry (summed shared
            cM within the population) -- recent-relatedness / endogamy structure.
  (b) SNP:  cluster individuals on pruned-SNP genotype PCA -- ancestry structure.

For each population and each method it runs PCA -> K-means (K chosen by silhouette,
K=1 if no structure) and keeps the LARGEST cluster.  Outputs two label files
(same 2-column '<sample> <POP>' format as the input) that SLICE the data when
passed to generate_stan_data.py --labels:

  <out>_ibd.txt   individuals selected by IBD clustering
  <out>_snp.txt   individuals selected by SNP clustering
  <out>_clusters.png   PC1-PC2 scatter per population/method (diagnostic)
  <out>_summary.json   cluster sizes and selections

Then regenerate + fit, e.g.:
  PY=/opt/miniconda3/envs/genetics_env/bin/python
  $PY generate_stan_data.py --labels <out>_snp.txt --pop-order AFR EAS EUR SAS \
      --ibd-method file --ibd-glob 'hapibd_merged/*.ibd.gz' \
      --vcf-glob '<pruned vcfs>' --genmap-dir merged_pruned/genmap \
      --bins-uniform 2 20.5 0.5 --min-maf 0.05 --out real_4pop_allchr_snpclust
  $PY fit_dense.py            # after pointing DATASETS at the sliced prefix
"""
import os
import glob
import gzip
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------
# IO
# ----------------------------------------------------------------------
def load_labels(path):
    """-> (ordered sample list, {sample: pop})."""
    samples, s2p = [], {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                samples.append(p[0]); s2p[p[0]] = p[1]
    return samples, s2p


def vcf_sample_order(vcf):
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.split()[9:]
    return []


def load_ibd_coancestry(ibd_glob, pop_members):
    """Per population, summed shared cM between each pair of its individuals.

    pop_members : {pop: [sample, ...]}   ->   {pop: n x n float matrix}
    """
    idx = {p: {s: i for i, s in enumerate(members)}
           for p, members in pop_members.items()}
    sample_pop = {s: p for p, members in pop_members.items() for s in members}
    M = {p: np.zeros((len(members), len(members)))
         for p, members in pop_members.items()}
    for f in sorted(glob.glob(ibd_glob)):
        with gzip.open(f, "rt") as fh:
            for line in fh:
                c = line.split()
                s1, s2, L = c[0], c[2], float(c[7])
                p = sample_pop.get(s1)
                if p is None or sample_pop.get(s2) != p:
                    continue                       # cross-pop or unlisted
                a, b = idx[p][s1], idx[p][s2]
                M[p][a, b] += L; M[p][b, a] += L
    return M


def load_genotypes(vcf_glob, want, max_snps, stride):
    """Dosage matrix for `want` samples across up-to max_snps pruned SNPs.

    Returns (kept_sample_list, D)  with D of shape (n_samples, n_snps).
    """
    vcfs = [v for v in sorted(glob.glob(vcf_glob)) if "chr9" not in os.path.basename(v)]
    order = vcf_sample_order(vcfs[0])
    cols = [i for i, s in enumerate(order) if s in want]
    kept = [order[i] for i in cols]
    rows = []
    seen = 0
    for v in vcfs:
        try:
            fh = gzip.open(v, "rt")
        except OSError:
            continue
        with fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                seen += 1
                if seen % stride:
                    continue
                gt = line.rstrip("\n").split("\t")[9:]
                d = np.empty(len(cols), dtype=np.float32)
                for k, ci in enumerate(cols):
                    s = gt[ci]
                    d[k] = np.nan if s[0] == "." else (s[0] == "1") + (s[2] == "1")
                rows.append(d)
                if len(rows) >= max_snps:
                    break
        if len(rows) >= max_snps:
            break
    D = np.asarray(rows, dtype=np.float32).T          # samples x snps
    # mean-impute missing per SNP
    colmean = np.nanmean(D, axis=0)
    inds = np.where(np.isnan(D))
    D[inds] = np.take(colmean, inds[1])
    return kept, D


# ----------------------------------------------------------------------
# PCA + clustering
# ----------------------------------------------------------------------
def pca(X, k):
    """Standardize columns, return top-k principal-component scores."""
    Xc = X - X.mean(0)
    sd = Xc.std(0); sd[sd == 0] = 1.0
    Xs = Xc / sd
    U, S, _ = np.linalg.svd(Xs, full_matrices=False)
    k = min(k, U.shape[1])
    return U[:, :k] * S[:k]


def silhouette(X, lab):
    """Mean silhouette over points (euclidean, in feature space)."""
    ks = np.unique(lab)
    if ks.size < 2:
        return -1.0
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    sil = np.zeros(len(X))
    for i in range(len(X)):
        same = lab == lab[i]; same[i] = False
        a = D[i, same].mean() if same.any() else 0.0
        b = min(D[i, lab == c].mean() for c in ks if c != lab[i])
        sil[i] = 0.0 if max(a, b) == 0 else (b - a) / max(a, b)
    return float(sil.mean())


def strip_outliers(feat, n_mad=6.0, ndim=5):
    """Robust inlier mask: drop points > n_mad robust-z on any top PC.

    Uses per-dimension median / MAD on the top `ndim` PCs, so a single extreme
    sample (e.g. a bad genotype vector) can't dominate the K-means / silhouette.
    """
    d = min(ndim, feat.shape[1])
    X = feat[:, :d]
    med = np.median(X, 0)
    mad = np.median(np.abs(X - med), 0) * 1.4826 + 1e-9
    z = np.abs(X - med) / mad
    return z.max(1) <= n_mad


def choose_clusters(feat, kmax, sil_thresh, seed, min_frac=0.10, n_mad=6.0):
    """Robust PCA-space K-means.

    (1) strip extreme outliers (label -1 = dropped);
    (2) K-means for K=2..kmax on the inliers, but only accept a split whose
        SMALLEST cluster is >= min_frac of the inliers (kills trivial
        outlier-only splits like 486/3);
    (3) pick the K with best silhouette; K=1 (all inliers one cluster) if none
        beats sil_thresh.

    Returns (labels over all n points with -1 for outliers, best_k, best_sil).
    """
    n = len(feat)
    labels = np.full(n, -1)
    inlier = strip_outliers(feat, n_mad)
    idx = np.where(inlier)[0]
    fin = feat[inlier]
    if len(fin) < 10:
        labels[inlier] = 0
        return labels, 1, 0.0
    min_size = max(5, int(min_frac * len(fin)))
    best = (None, 1, -1.0)
    for K in range(2, min(kmax, len(fin) // min_size) + 1):
        _, lab = kmeans2(fin, K, minit="++", seed=seed, missing="warn")
        if len(np.unique(lab)) < K:                       # collapsed
            continue
        if np.bincount(lab, minlength=K).min() < min_size:  # trivial split
            continue
        s = silhouette(fin, lab)
        if s > best[2]:
            best = (lab, K, s)
    if best[0] is None or best[2] < sil_thresh:           # no real structure
        labels[inlier] = 0
        return labels, 1, (best[2] if best[0] is not None else 0.0)
    labels[idx] = best[0]
    return labels, best[1], best[2]


def largest_cluster_mask(lab):
    """Largest real cluster (ignores outliers labelled -1)."""
    valid = lab[lab >= 0]
    if valid.size == 0:
        return np.zeros(len(lab), bool)
    vals, cnts = np.unique(valid, return_counts=True)
    return lab == vals[np.argmax(cnts)]


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", default=os.path.join(_THIS_DIR, "population_labels_4pop.txt"))
    ap.add_argument("--pop-order", nargs="+", default=["AFR", "EAS", "EUR", "SAS"])
    ap.add_argument("--ibd-glob", default=os.path.join(_THIS_DIR, "hapibd_merged", "*.ibd.gz"))
    ap.add_argument("--vcf-glob", default=os.path.join(_THIS_DIR, "merged_pruned", "vcf", "*.vcf.gz"))
    ap.add_argument("--out", default=os.path.join(_THIS_DIR, "cluster4pop"))
    ap.add_argument("--n-pc", type=int, default=10)
    ap.add_argument("--kmax", type=int, default=5)
    ap.add_argument("--sil-thresh", type=float, default=0.10,
                    help="Min mean silhouette to accept >1 cluster (else keep all).")
    ap.add_argument("--min-frac", type=float, default=0.10,
                    help="Min cluster size as a fraction of inliers (rejects "
                         "trivial outlier-only splits).")
    ap.add_argument("--n-mad", type=float, default=6.0,
                    help="Robust-z (MAD) threshold for dropping PC outliers.")
    ap.add_argument("--max-snps", type=int, default=20000)
    ap.add_argument("--snp-stride", type=int, default=10,
                    help="Take every Nth pruned SNP (spread across the genome).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    np.random.seed(args.seed)

    samples, s2p = load_labels(args.labels)
    pops = args.pop_order
    pop_members = {p: [s for s in samples if s2p.get(s) == p] for p in pops}
    for p in pops:
        print(f"[load] {p}: {len(pop_members[p])} individuals")

    # ---- features ----
    print("[ibd] building co-ancestry matrices ...")
    coanc = load_ibd_coancestry(args.ibd_glob, pop_members)
    print(f"[snp] loading genotypes (<= {args.max_snps} SNPs, stride {args.snp_stride}) ...")
    want = set(samples)
    gsamp, D = load_genotypes(args.vcf_glob, want, args.max_snps, args.snp_stride)
    gidx = {s: i for i, s in enumerate(gsamp)}
    print(f"[snp] genotype matrix {D.shape} (samples x snps)")

    methods = {}
    summary = {}
    for method, feat_fn in [
        ("ibd", lambda p: pca(np.log1p(coanc[p]), args.n_pc)),
        ("snp", lambda p: pca(D[[gidx[s] for s in pop_members[p] if s in gidx]], args.n_pc)),
    ]:
        selected = {}          # pop -> list of kept samples
        scores = {}            # pop -> (feat, labels, keepmask, members)
        for p in pops:
            members = ([s for s in pop_members[p] if s in gidx]
                       if method == "snp" else pop_members[p])
            feat = feat_fn(p)
            lab, K, sil = choose_clusters(feat, args.kmax, args.sil_thresh, args.seed,
                                          min_frac=args.min_frac, n_mad=args.n_mad)
            keep = largest_cluster_mask(lab)
            sel = [members[i] for i in range(len(members)) if keep[i]]
            selected[p] = sel
            scores[p] = (feat, lab, keep, members)
            uv, uc = np.unique(lab, return_counts=True)
            n_out = int(uc[uv == -1].sum()) if (uv == -1).any() else 0
            sizes = [int(c) for v, c in zip(uv, uc) if v >= 0]
            summary.setdefault(p, {})[method] = {
                "n_total": len(members), "K": int(K), "silhouette": round(sil, 3),
                "cluster_sizes": sizes, "n_outliers": n_out, "n_selected": len(sel)}
            print(f"[{method}] {p}: n={len(members)} K={K} sil={sil:.2f} "
                  f"sizes={sizes} outliers={n_out} -> keep {len(sel)}")
        methods[method] = (selected, scores)

    # ---- write sliced label files ----
    for method, (selected, _) in methods.items():
        out = f"{args.out}_{method}.txt"
        with open(out, "w") as f:
            for p in pops:
                for s in selected[p]:
                    f.write(f"{s} {p}\n")
        n = sum(len(selected[p]) for p in pops)
        print(f"[save] {out}  ({n} individuals)")

    # ---- diagnostic PCA scatter ----
    fig, axes = plt.subplots(2, len(pops), figsize=(4 * len(pops), 8), squeeze=False)
    for r, method in enumerate(["ibd", "snp"]):
        _, scores = methods[method]
        for c, p in enumerate(pops):
            ax = axes[r][c]
            feat, lab, keep, _ = scores[p]
            ax.scatter(feat[~keep, 0], feat[~keep, 1], s=8, c="#bbbbbb", label="dropped")
            ax.scatter(feat[keep, 0], feat[keep, 1], s=8, c="#c0392b", label="kept")
            ax.set_title(f"{method.upper()}  {p}  (keep {keep.sum()}/{len(keep)})", fontsize=10)
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    axes[0][0].legend(fontsize=8, loc="best")
    fig.suptitle("Within-population clustering: kept (red) vs dropped (grey)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{args.out}_clusters.png", dpi=150, bbox_inches="tight")
    print(f"[save] {args.out}_clusters.png")

    json.dump(summary, open(f"{args.out}_summary.json", "w"), indent=2)
    print(f"[save] {args.out}_summary.json")


if __name__ == "__main__":
    main()
