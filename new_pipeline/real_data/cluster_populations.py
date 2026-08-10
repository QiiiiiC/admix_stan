"""Find a homogeneous sub-cluster within each superpopulation, two ways.

The coalescent model assumes any two individuals *within* a leaf population can
coalesce (panmixia).  That is violated when a 1000G superpopulation pools several
subpopulations (e.g. AFR = YRI/LWK/GWD/MSL/ESN + admixed ASW/ACB).  This script
detects that substructure and picks the dominant homogeneous cluster per
population, by two independent signals:

  (a) IBD:  cluster individuals on their pairwise IBD co-ancestry (summed shared
            cM within the population) -- recent-relatedness / endogamy structure.
  (b) SNP:  cluster individuals on pruned-SNP genotype PCA -- ancestry structure.

The two signals need different machinery, because they are different objects:

  SNP:  genotypes are a real individuals x markers DATA matrix, so the standard
        popgen recipe applies -- standardize markers, PCA (Patterson, Price &
        Reich 2006; scores kept at natural scale, NOT whitened), then K-means
        with K chosen by silhouette.
  IBD:  co-ancestry is an n x n SIMILARITY matrix with no coordinate space, so
        it is clustered as a GRAPH: Louvain modularity, normalized-Laplacian
        spectral clustering, and (for reference) DBSCAN on PC scores.  All three
        run every time; --ibd-cluster picks which writes the labels.

Which cluster then REPRESENTS the population is a two-step rule (--select
cohesive, the default): discard every cluster smaller than --min-cluster-size,
and among the survivors take the most internally cohesive one.  The size floor
comes first because cohesion is maximised by tiny cliques -- a pair of cousins
is the tightest "population" in any dataset, and useless downstream.  Taking the
largest cluster instead (--select largest) is the old behaviour; it is the wrong
default here, because when a superpopulation splits the biggest piece is often
the loose admixed remainder rather than the homogeneous core.

Cohesion is measured against each signal's own null, since the two objects are
different (see below):

  IBD:  connectivity -- the fraction of within-cluster pairs sharing any IBD.
        A panmictic leaf is one where every pair can coalesce recently, so this
        measures the model assumption directly.  Note it is NOT mean co-ancestry
        enrichment; that was tried and is anti-correlated with fit quality, for
        the reason given on ibd_cohesion().
  SNP:  mean within-cluster silhouette in PC space -- the same criterion that
        chose K, restricted to one cluster.

If the 1000G panel file is present, every partition is scored against the
reported subpopulations (ARI/AMI) -- see the [truth] block, and note that the
SNP route currently beats all IBD routes.

  <out>_ibd.txt   individuals selected by IBD clustering
  <out>_snp.txt   individuals selected by SNP clustering
  <out>_clusters.png   IBD: co-ancestry enrichment heatmaps reordered by cluster
                       (red diagonal blocks = real communities); SNP: PC scatter
  <out>_summary.json   cluster sizes, selections, agreement with the panel

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
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
import networkx as nx
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics import (silhouette_score, adjusted_rand_score,
                             adjusted_mutual_info_score)
from sklearn.neighbors import NearestNeighbors

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


def base_id(s):
    """Strip the plink FID_IID doubling: 'HG00096_HG00096' -> 'HG00096'.

    The IBD files carry plain sample IDs while the pruned VCFs (written through
    plink) carry the doubled form, so every join between the two sources has to
    go through one canonical key.  1000G IDs contain no underscore themselves.
    """
    return s.split("_")[0]


def vcf_sample_order(vcf):
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.split()[9:]
    return []


def load_ibd_coancestry(ibd_glob, pop_members, min_cm=2.0, max_cm=np.inf,
                        genome_cm=None):
    """Per population, the pairwise IBD SHARING FRACTION between its individuals.

    W_ij = (total cM shared by i and j in [min_cm, max_cm)) / genome_cm, i.e. the
    fraction of the autosome the pair shares IBD.  With genome_cm=None the raw
    summed cM is returned instead.

    The normalisation is cosmetic for the clustering itself -- Louvain's
    modularity, the normalized Laplacian and the connectivity score are all
    invariant to a global rescaling of W, so dividing by a constant cannot move
    a single cluster boundary.  It is worth doing anyway because it puts the
    matrix in units you can read: W_ij = 0.01 means "these two share 1% of their
    genome", which is roughly a second cousin.

    Only segments with min_cm <= length < max_cm contribute.  Segment length is a
    clock: a segment of L cM has expected TMRCA ~ 100/(2L) generations, so the
    bin selects the time depth the co-ancestry matrix describes.  >10 cM is
    pedigree relatedness; the SHORT end is where subpopulation structure lives.

    The floor defaults to 2 cM.  An earlier version of this script used 1 cM on
    the OLD pipeline (reference-free Beagle phasing, full sequence density),
    where it scored better against the reported subpopulations -- but on that
    data the sub-2 cM class was largely phase-switch debris, so the floor was
    compensating for a defect that the high_cov pipeline has since removed.  On
    the current data the two floors are compared directly in the run log; 2 cM
    is the honest default because it is the length class HapNe itself fits and
    the one whose detection is reliable.

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
                if L < min_cm or L >= max_cm:
                    continue
                p = sample_pop.get(s1)
                if p is None or sample_pop.get(s2) != p:
                    continue                       # cross-pop or unlisted
                a, b = idx[p][s1], idx[p][s2]
                M[p][a, b] += L; M[p][b, a] += L
    if genome_cm:
        for p in M:
            M[p] /= genome_cm
    return M


def load_genotypes(vcf_glob, want, max_snps, stride):
    """Dosage matrix for `want` samples across up-to max_snps pruned SNPs.

    `want` holds canonical (base) IDs; VCF columns are matched through base_id
    and returned under the canonical name, so callers never see the FID_IID form.

    Returns (kept_sample_list, D)  with D of shape (n_samples, n_snps).
    """
    vcfs = [v for v in sorted(glob.glob(vcf_glob)) if "chr9" not in os.path.basename(v)]
    order = vcf_sample_order(vcfs[0])
    cols = [i for i, s in enumerate(order) if base_id(s) in want]
    kept = [base_id(order[i]) for i in cols]
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
    """Normalize FEATURES, return top-k PC scores at their natural scale.

    This is the standard population-genetics PCA (Patterson, Price & Reich 2006;
    smartpca / plink --pca): each marker/column is centred and divided by its own
    standard deviation, then the scores are U*S -- eigenvector times sqrt of
    eigenvalue.

    Do NOT rescale the PC *scores* to unit variance.  The relative scale of the
    axes is the signal: an axis' singular value is how much population
    differentiation it carries, so it should down-weight itself as it decays into
    the noise bulk.  Whitening throws that away and promotes every retained noise
    axis to the same weight as PC1 -- measured here, it destroyed all structure
    by k=5, while the unwhitened scores stay stable out to k=20.

    Returns (n x k scores, explained-variance ratio of every PC).
    """
    Xc = X - X.mean(0)
    sd = Xc.std(0); sd[sd == 0] = 1.0
    Xs = Xc / sd
    U, S, _ = np.linalg.svd(Xs, full_matrices=False)
    ev = (S ** 2) / (S ** 2).sum()
    k = min(k, U.shape[1])
    return U[:, :k] * S[:k], ev


def strip_outliers(feat, n_mad=6.0, ndim=5):
    """Robust inlier mask: drop points > n_mad robust-z on any of the top PCs.

    K-means has no noise class, so one extreme sample (a bad genotype vector, a
    lone recent migrant) drags a centroid onto itself.  Per-dimension median/MAD
    on the leading PCs removes those before fitting.
    """
    d = min(ndim, feat.shape[1])
    X = feat[:, :d]
    med = np.median(X, 0)
    mad = np.median(np.abs(X - med), 0) * 1.4826 + 1e-9
    return (np.abs(X - med) / mad).max(1) <= n_mad


def kmeans_pcs(feat, kmax=8, sil_thresh=0.10, seed=0, min_frac=0.10, n_mad=6.0):
    """PCA-space K-means with K chosen by silhouette -- the standard SNP recipe.

    Genotypes are a real data matrix, so the pipeline is the textbook one:
    standardize markers -> PCA -> K-means on the unwhitened PC scores.

    Two guards, both needed because K-means partitions EVERY point:
      (1) strip MAD outliers first (labelled -1), so they cannot capture centroids;
      (2) reject any K whose smallest cluster is under min_frac of the inliers.
    Without (2) silhouette runs away to large K, because a cluster holding one
    isolated outlier scores near 1.0 and drags the mean up -- unguarded this gave
    AFR K=7 with clusters of size 4 and 2, and shrank the kept set to 231/661.

    Returns (labels with -1 for outliers, K, silhouette).
    """
    n = len(feat)
    labels = np.full(n, -1)
    inlier = strip_outliers(feat, n_mad)
    idx = np.where(inlier)[0]
    fin = feat[inlier]
    if len(fin) < 20:
        labels[inlier] = 0
        return labels, 1, 0.0
    min_size = max(5, int(min_frac * len(fin)))
    best = (None, 1, -1.0)
    for K in range(2, min(kmax, len(fin) // min_size) + 1):
        lab = KMeans(n_clusters=K, n_init=10, random_state=seed).fit_predict(fin)
        if len(np.unique(lab)) < K:                        # collapsed
            continue
        if np.bincount(lab, minlength=K).min() < min_size:  # trivial split
            continue
        s = float(silhouette_score(fin, lab))
        if s > best[2]:
            best = (lab, K, s)
    if best[0] is None or best[2] < sil_thresh:
        labels[inlier] = 0
        return labels, 1, (best[2] if best[0] is not None else 0.0)
    order = np.argsort(-np.bincount(best[0], minlength=best[1]))
    remap = np.empty(best[1], int); remap[order] = np.arange(best[1])
    labels[idx] = remap[best[0]]
    return labels, best[1], best[2]


def knee_eps(feat, min_samples):
    """DBSCAN eps from the k-distance knee (Ester et al. 1996).

    Sort every point's distance to its `min_samples`-th nearest neighbour; the
    curve is flat over the dense core and turns up sharply where the sparse tail
    begins.  Take the knee = the point furthest from the chord joining the two
    ends of the sorted curve.
    """
    k = min(min_samples, len(feat) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(feat)   # +1: self is neighbour 0
    d = np.sort(nn.kneighbors(feat)[0][:, k])
    x = np.arange(len(d), dtype=float)
    # distance from each (x, d) to the chord (x0,d0)-(x1,d1), unnormalised
    dx, dy = x[-1] - x[0], d[-1] - d[0]
    norm = np.hypot(dx, dy)
    if norm < 1e-12:
        return float(d[-1]) if d[-1] > 0 else 1.0
    off = np.abs(dy * (x - x[0]) - dx * (d - d[0])) / norm
    eps = float(d[int(np.argmax(off))])
    return eps if eps > 0 else float(d[d > 0].min() if (d > 0).any() else 1.0)


def cluster_dbscan(feat, min_samples, eps=None):
    """DBSCAN on an embedding (Laplacian for IBD, PC scores for SNP).

    Unlike K-means there is no K to choose and no separate outlier pass: points
    in low-density regions come back as noise (-1), which is exactly the
    admixed / related stragglers we want out of the panmictic leaf.

    For the IBD route `feat` is the normalized-Laplacian embedding, not a PCA of
    the co-ancestry matrix.  Running PCA on W treats each ROW of a similarity
    matrix as if it were a feature vector for that individual, which is not what
    it is; the Laplacian embedding is the coordinate space the graph actually
    induces.

    Returns (labels with -1 = noise, n_clusters, mean silhouette over non-noise,
    eps actually used).
    """
    n = len(feat)
    if n < min_samples + 1:
        return np.zeros(n, int), 1, 0.0, float("nan")
    e = knee_eps(feat, min_samples) if eps is None else eps
    lab = DBSCAN(eps=e, min_samples=min_samples).fit_predict(feat)
    ks = np.unique(lab[lab >= 0])
    if ks.size == 0:                                  # eps too tight: all noise
        return np.zeros(n, int), 1, 0.0, e
    sil = 0.0
    if ks.size > 1:
        m = lab >= 0
        sil = float(silhouette_score(feat[m], lab[m]))
    return lab, int(ks.size), sil, e


# ----------------------------------------------------------------------
# IBD graph clustering (operates on the co-ancestry matrix directly)
# ----------------------------------------------------------------------
def louvain_clusters(W, resolution=1.0, seed=0):
    """Louvain modularity communities on the weighted IBD graph.

    Nodes are individuals, edge weight = total shared cM.  Modularity is

        Q = 1/2m * sum_ij [ w_ij - k_i k_j / 2m ] delta(c_i, c_j)

    and the k_i k_j / 2m term is a configuration-model null: it asks whether i
    and j share MORE than expected given how much each shares overall.  That
    degree correction is why the raw matrix needs no hand-scaling -- individuals
    differ in total sharing through phasing quality and through genuine endogamy,
    and row-normalising would delete the second (real) effect along with the
    first.  Han et al. 2017 (Nat Commun) cluster 770k genomes this way.

    Returns integer labels (no noise class; every node joins a community).
    """
    G = nx.Graph()
    G.add_nodes_from(range(len(W)))
    iu = np.triu_indices(len(W), 1)
    w = W[iu]
    nz = w > 0
    G.add_weighted_edges_from(zip(iu[0][nz], iu[1][nz], w[nz]))
    comms = nx.community.louvain_communities(G, weight="weight",
                                             resolution=resolution, seed=seed)
    comms = sorted(comms, key=len, reverse=True)
    lab = np.zeros(len(W), int)
    for c, nodes in enumerate(comms):
        for i in nodes:
            lab[i] = c
    return lab


def knn_sparsify(W, n_neighbors=20):
    """Keep each node's strongest n_neighbors edges, then symmetrise (OR-rule).

    Spectral clustering needs a SPARSE affinity.  On a fully connected graph the
    normalized Laplacian has one zero eigenvalue and the rest spread smoothly, so
    the spectrum shows no block structure and the eigengap is uninformative --
    every population came back as a single cluster before this step.  von Luxburg
    (2007) recommends a k-NN graph for exactly this reason.
    """
    n = len(W)
    k = min(n_neighbors, n - 1)
    keep = np.zeros_like(W, dtype=bool)
    nb = np.argpartition(-W, k, axis=1)[:, :k]
    keep[np.arange(n)[:, None], nb] = True
    keep |= keep.T                          # OR-rule symmetrisation
    return np.where(keep, W, 0.0)


def laplacian_embedding(W, k, n_neighbors=20):
    """Row-normalised bottom-k eigenvectors of L_sym = I - D^-1/2 W D^-1/2.

    This is the coordinate space a similarity matrix actually has.  The raw
    co-ancestry matrix has none -- W_ij is one number per pair, not a position
    -- so anything that needs distances (DBSCAN, k-means) has to be given an
    embedding first, and the normalized Laplacian is the principled one: its
    bottom eigenvectors are the relaxed solution to normalized min-cut (Shi &
    Malik 2000; von Luxburg 2007), and the D^-1/2 on both sides applies the same
    degree correction that Louvain's null model applies.

    Used by both spectral_clusters (k-means on it) and cluster_dbscan (density
    clustering on it).
    """
    n = len(W)
    W = knn_sparsify(W, n_neighbors)
    d = W.sum(1)
    d[d <= 0] = 1e-12
    dis = 1.0 / np.sqrt(d)
    L = np.eye(n) - (W * dis[:, None]) * dis[None, :]
    vals, vecs = np.linalg.eigh(L)
    k = max(2, min(k, n - 1))
    U = vecs[:, :k]
    rn = np.linalg.norm(U, axis=1, keepdims=True); rn[rn < 1e-12] = 1.0
    return U / rn, vals


def spectral_clusters(W, kmax=8, n_neighbors=20, seed=0):
    """Normalized-Laplacian spectral clustering (Ng-Jordan-Weiss 2002).

    L_sym = I - D^-1/2 W D^-1/2 on a k-NN-sparsified affinity; the number of
    clusters comes from the eigengap, then k-means on the row-normalised bottom-k
    eigenvectors.  D^-1/2 on both sides is the same degree correction modularity
    gets from its null model.

    This is the correct eigen-method for a similarity matrix -- note it is NOT
    what an SVD of the raw co-ancestry matrix computes, since that treats rows of
    a similarity matrix as if they were feature vectors.

    Returns (labels, k chosen by the eigengap).
    """
    n = len(W)
    _, vals = laplacian_embedding(W, 2, n_neighbors)
    # Eigengap over k >= 2.  The k=1 gap (lambda_1 - lambda_0) is always large --
    # lambda_0 is the trivial 0 of a connected graph -- so including it would
    # answer "one cluster" every time regardless of the structure present.
    top = min(kmax + 1, n - 1)
    gaps = vals[2:top] - vals[1:top - 1]
    if gaps.size == 0:
        return np.zeros(n, int), 1
    k = int(np.argmax(gaps)) + 2
    U, _ = laplacian_embedding(W, k, n_neighbors)
    lab = KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(U)
    # relabel so cluster 0 is the largest, matching louvain_clusters
    order = np.argsort(-np.bincount(lab, minlength=k))
    remap = np.empty(k, int); remap[order] = np.arange(k)
    return remap[lab], k


def plot_coancestry_blocks(ax, W, lab, title, chosen=None):
    """Reordered co-ancestry heatmap -- the honest view of an n x n matrix.

    The IBD data has no coordinate space: it is one similarity per pair, so a
    2-D scatter of PC scores is a projection of something that was never a point
    cloud.  Sorting rows/columns by cluster and showing the matrix directly is
    what fineSTRUCTURE plots.

    What is plotted is ENRICHMENT, log2 of W_ij / (k_i k_j / 2m), not raw sharing.
    Raw co-ancestry has almost no visible contrast: every pair in a population
    shares a large common baseline, and individuals differ several-fold in total
    sharing, so rows stripe by degree and the blocks wash out.  Dividing by the
    configuration-model expectation -- the same null Louvain's modularity uses --
    removes both effects and leaves "more/less than expected for these two".  Red
    diagonal blocks against a blue off-diagonal is a real partition.
    """
    order = np.concatenate([np.where(lab == c)[0] for c in
                            sorted(set(lab), key=lambda c: -(lab == c).sum())])
    k = W.sum(1); m = W.sum() / 2.0
    E = np.outer(k, k) / (2 * m) if m > 0 else np.ones_like(W)
    E[E <= 0] = 1e-12
    R = np.log2(np.clip(W / E, 1e-3, None))[np.ix_(order, order)]
    np.fill_diagonal(R, 0.0)                    # self-sharing is not data
    # Pair-level values are very noisy, and at ~500x500 in a small panel that
    # speckle hides the blocks.  Average into tiles so each pixel is a mean over
    # many pairs -- the block signal survives, the pair noise does not.
    n = len(R); tgt = 120
    if n > tgt:
        t = int(np.ceil(n / tgt)); pad = (-n) % t
        Rp = np.pad(R, ((0, pad), (0, pad)), constant_values=np.nan)
        R = np.nanmean(Rp.reshape(Rp.shape[0] // t, t, -1, t), axis=(1, 3))
        scale = 1.0 / t
    else:
        scale = 1.0
    # Asymmetric limits centred on 0, NOT +/- one symmetric value.  The
    # enrichment distribution is strongly one-sided: most pairs in a population
    # share nothing, so their log2 ratio sits at the clip floor (~-10) while a
    # real diagonal block only reaches +1 or so.  A symmetric scale set from
    # |R| is therefore fixed by the empty pairs, the whole red half goes unused,
    # and every genuine community renders as near-white -- the blocks are there
    # but invisible.  Taking the 2nd/98th percentiles and pinning white to 0
    # (TwoSlopeNorm) spends each half of the colormap on the range that half
    # actually covers, so "more than expected" finally reads as red.
    lo = min(float(np.nanpercentile(R, 2)), -1e-3)
    hi = max(float(np.nanpercentile(R, 98)), 1e-3)
    norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi)
    im = ax.imshow(R, cmap="RdBu_r", norm=norm, interpolation="nearest")
    b = 0
    for c in sorted(set(lab), key=lambda c: -(lab == c).sum()):
        n_c = int((lab == c).sum())
        if chosen is not None and c == chosen:
            # green box = the cluster the selection rule kept for this algorithm
            blo, bhi = b * scale - .5, (b + n_c) * scale - .5
            ax.add_patch(Rectangle((blo, blo), bhi - blo, bhi - blo, fill=False,
                                   edgecolor="#00b050", lw=2.5, zorder=5))
        b += n_c
        if b < len(lab):
            ax.axhline(b * scale - .5, color="#111111", lw=1.0)
            ax.axvline(b * scale - .5, color="#111111", lw=1.0)
    ax.set_title(title, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.02,
                            ticks=[lo, 0, hi])
    cb.ax.set_yticklabels(["less\nthan expected", "as\nexpected", "more\nthan expected"],
                          fontsize=6)
    cb.ax.tick_params(length=0)


def largest_cluster_mask(lab):
    """Largest real cluster (ignores outliers labelled -1)."""
    valid = lab[lab >= 0]
    if valid.size == 0:
        return np.zeros(len(lab), bool)
    vals, cnts = np.unique(valid, return_counts=True)
    return lab == vals[np.argmax(cnts)]


# ----------------------------------------------------------------------
# "how close is this cluster" -- the score the selection rule maximises
# ----------------------------------------------------------------------
def ibd_cohesion(W, lab):
    """Per-cluster CONNECTIVITY: fraction of within-cluster pairs sharing any IBD.

        cohesion(C) = |{ i<j in C : W_ij > 0 }| / |{ i<j in C }|

    This is the quantity the downstream model actually assumes.  A coalescent
    leaf is panmictic: any two members can find a recent common ancestor, so
    every pair should show some co-ancestry.  Connectivity measures exactly that
    and nothing else.

    It replaces the obvious alternative -- mean co-ancestry enrichment,
    mean(W_ij) / mean(k_i k_j / 2m) -- which was measured here and is ANTI-
    correlated with whether the resulting sample fits.  IBD sharing per pair is
    heavy-tailed, so a mean is set by a handful of cryptic relatives: it rewards
    a cluster containing a few very-close pairs over one where everybody is
    modestly related, which is the opposite of exchangeability.  Concretely, on
    these four superpopulations mean enrichment ranked FIN LAST inside EUR
    (1.56, below TSI's 7.43) -- and FIN is the single population that reproduces
    the published HapNe-IBD curve.  Connectivity ranks it first at 98.7%, and
    picks LWK in AFR and GIH in SAS, the three fits that land near the ~1e4
    ancestral Ne while every low-connectivity pick rails to 1e5-1e6.

    Combined with the --min-cluster-size floor this is not gameable by small
    cliques: a pair of siblings scores 100% but never clears the floor.
    """
    out = {}
    for c in sorted({int(v) for v in lab if v >= 0}):
        idx = np.where(lab == c)[0]
        if len(idx) < 2:
            out[c] = 0.0
            continue
        iu = np.triu_indices(len(idx), 1)
        out[c] = float((W[np.ix_(idx, idx)][iu] > 0).mean())
    return out


def snp_cohesion(feat, lab):
    """Per-cluster mean silhouette in PC space.

    Silhouette is already the criterion that chose K, so restricting it to one
    cluster keeps the SNP route internally consistent: s(i) contrasts i's mean
    distance to its own cluster against its mean distance to the nearest other
    one, so a high-scoring cluster is both tight AND well separated -- exactly
    "the most distinct homogeneous group in this population".

    Undefined for a single cluster (there is no other cluster to contrast
    against); that case scores 1.0, since it is selected regardless.
    """
    ok = lab >= 0
    ks = {int(v) for v in lab[ok]}
    if len(ks) < 2:
        return {c: 1.0 for c in ks}
    from sklearn.metrics import silhouette_samples
    s = silhouette_samples(feat[ok], lab[ok])
    sub_ = lab[ok]
    return {c: float(s[sub_ == c].mean()) for c in ks}


def select_cluster(lab, cohesion, min_size, mode="cohesive"):
    """Which cluster represents the population.

    'cohesive' (default): drop clusters under min_size, then take the highest
    cohesion among what is left.  The size floor has to come first -- cohesion
    is monotonically easier to achieve as a cluster shrinks, and unbounded it
    would always return the tightest handful of close relatives, which is both
    unrepresentative and too small for anything downstream to fit.

    'largest': the biggest cluster, ignoring cohesion.

    Falls back to the largest cluster when nothing clears min_size, so a
    population is never left with an empty selection.

    Returns (boolean mask, chosen label, reason string).
    """
    valid = lab[lab >= 0]
    if valid.size == 0:
        return np.zeros(len(lab), bool), None, "no clusters"
    vals, cnts = np.unique(valid, return_counts=True)
    if mode == "largest":
        c = int(vals[np.argmax(cnts)])
        return lab == c, c, "largest"
    ok = [(int(v), int(n)) for v, n in zip(vals, cnts) if n >= min_size]
    if not ok:
        c = int(vals[np.argmax(cnts)])
        return lab == c, c, f"fallback: no cluster >= {min_size}"
    c = max(ok, key=lambda t: cohesion.get(t[0], -np.inf))[0]
    return lab == c, c, f"most cohesive of {len(ok)} clusters >= {min_size}"


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", default=os.path.join(_THIS_DIR, "population_labels_4pop.txt"))
    ap.add_argument("--pop-order", nargs="+", default=["AFR", "EAS", "EUR", "SAS"])
    ap.add_argument("--ibd-glob", default=os.path.join(_THIS_DIR, "hapibd_merged", "*.ibd.gz"))
    ap.add_argument("--vcf-glob", default=os.path.join(_THIS_DIR, "merged_pruned", "vcf", "*.vcf.gz"))
    ap.add_argument("--out", default=os.path.join(_THIS_DIR, "cluster4pop"))
    ap.add_argument("--n-pc", type=int, default=4,
                    help="PCs kept, at their natural (unwhitened) scale. Scores "
                         "self-weight by singular value, so this is not critical; "
                         "the per-population variance spectrum is printed so you "
                         "can see where the noise bulk starts.")
    ap.add_argument("--kmax", type=int, default=8,
                    help="Max K for the SNP K-means sweep (best silhouette wins).")
    ap.add_argument("--sil-thresh", type=float, default=0.10,
                    help="Min mean silhouette to accept K>1 for the SNP K-means; "
                         "below it the population is called homogeneous.")
    ap.add_argument("--min-samples", type=int, default=10,
                    help="DBSCAN min_samples: neighbourhood size a point needs to "
                         "be a core point. Doubles as the smallest cluster size.")
    ap.add_argument("--eps", type=float, default=None,
                    help="DBSCAN eps in PC-score units. Default: per-population "
                         "k-distance knee.")
    ap.add_argument("--ibd-min-cm", type=float, default=2.0,
                    help="Only IBD segments >= this length build the co-ancestry "
                         "graph (~TMRCA 100/2L generations; 2 cM ~ 25 gen). This "
                         "is the length class HapNe fits and the one hap-ibd "
                         "detects reliably on the high_cov data.")
    ap.add_argument("--genome-cm", type=float, default=3372.18,
                    help="Total autosomal genetic length used to turn summed "
                         "shared cM into a SHARING FRACTION. Default is the sum "
                         "of the 39 GRCh38 masked chromosome arms. Clustering is "
                         "invariant to this constant; it only sets the units.")
    ap.add_argument("--ibd-max-cm", type=float, default=np.inf,
                    help="Upper end of the segment-length bin (default open).")
    ap.add_argument("--panel", default=os.path.join(
                        _THIS_DIR, "integrated_call_samples_v3.20130502.ALL.panel"),
                    help="1000G panel file; if present, every partition is scored "
                         "against the reported subpopulation labels (ARI/AMI).")
    ap.add_argument("--ibd-cluster", choices=["louvain", "spectral", "dbscan"],
                    default="spectral",
                    help="Which IBD partition writes the label file. All three are "
                         "always computed and compared.")
    ap.add_argument("--resolution", type=float, default=1.0,
                    help="Louvain resolution; >1 gives more, smaller communities.")
    ap.add_argument("--knn", type=int, default=20,
                    help="Neighbours kept per node when sparsifying the IBD graph "
                         "for spectral clustering.")
    ap.add_argument("--select", choices=["cohesive", "largest"], default="cohesive",
                    help="Which cluster represents the population. 'cohesive': "
                         "among clusters of at least --min-cluster-size members, "
                         "the most internally cohesive (IBD: fraction of pairs "
                         "sharing any IBD; SNP: mean silhouette). 'largest': "
                         "biggest cluster.")
    ap.add_argument("--min-cluster-size", type=int, default=50,
                    help="Size floor applied before the cohesion comparison. "
                         "Needed because cohesion rises as a cluster shrinks, so "
                         "unbounded it returns a handful of close relatives. The "
                         "floor is also a downstream requirement: HapNe-IBD needs "
                         "enough within-population pairs to see long segments "
                         "(~100 individuals reproduced the published FIN curve).")
    ap.add_argument("--max-snps", type=int, default=20000)
    ap.add_argument("--snp-stride", type=int, default=10,
                    help="Take every Nth pruned SNP (spread across the genome).")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    np.random.seed(args.seed)

    samples, s2p = load_labels(args.labels)
    pops = args.pop_order
    truth = {}
    if args.panel and os.path.exists(args.panel):
        with open(args.panel) as f:
            next(f, None)
            for line in f:
                c = line.split()
                if len(c) >= 3:
                    truth[c[0]] = c[1]
        print(f"[panel] {len(truth)} reported subpopulation labels from {args.panel}")
    pop_members = {p: [s for s in samples if s2p.get(s) == p] for p in pops}
    for p in pops:
        print(f"[load] {p}: {len(pop_members[p])} individuals")

    # ---- features ----
    hi = "inf" if not np.isfinite(args.ibd_max_cm) else f"{args.ibd_max_cm:g}"
    print(f"[ibd] building co-ancestry matrices, segments [{args.ibd_min_cm:g}, {hi}) cM ...")
    coanc = load_ibd_coancestry(args.ibd_glob, pop_members,
                                args.ibd_min_cm, args.ibd_max_cm, args.genome_cm)
    print(f"[snp] loading genotypes (<= {args.max_snps} SNPs, stride {args.snp_stride}) ...")
    want = set(samples)
    gsamp, D = load_genotypes(args.vcf_glob, want, args.max_snps, args.snp_stride)
    gidx = {s: i for i, s in enumerate(gsamp)}
    print(f"[snp] genotype matrix {D.shape} (samples x snps)")

    def describe(lab):
        uv, uc = np.unique(lab, return_counts=True)
        n_out = int(uc[uv == -1].sum()) if (uv == -1).any() else 0
        sizes = sorted((int(c) for v, c in zip(uv, uc) if v >= 0), reverse=True)
        return sizes, n_out

    def compo(lab, coh, chosen):
        """'size(cohesion)' per cluster, marking the selected one with *."""
        uv, uc = np.unique(lab[lab >= 0], return_counts=True)
        o = np.argsort(-uc)
        return " ".join(
            f"{'*' if int(uv[i]) == chosen else ''}{uc[i]}"
            f"({coh.get(int(uv[i]), float('nan')):.2f})" for i in o)

    methods = {}
    summary = {}
    ibd_alt = {}               # pop -> {algo: labels}, for the comparison table
    ibd_pick = {}              # pop -> {algo: chosen label}, for the plot marks
    ibd_sel = {}               # algo -> {pop: [selected samples]}, one file each

    # ---- IBD: graph clustering on the co-ancestry matrix ----
    selected, scores = {}, {}
    for p in pops:
        members = pop_members[p]
        W = coanc[p].copy()
        np.fill_diagonal(W, 0.0)
        # One embedding, used by both eigen-methods: DBSCAN needs distances and
        # the Laplacian is where the graph's distances live (see
        # laplacian_embedding).  n_pc components is plenty -- the structure sits
        # in the bottom few eigenvectors.
        emb, _ = laplacian_embedding(W, args.n_pc, args.knn)
        parts = {
            "louvain": louvain_clusters(W, args.resolution, args.seed),
            "spectral": spectral_clusters(W, n_neighbors=args.knn, seed=args.seed)[0],
            "dbscan": cluster_dbscan(emb, args.min_samples, args.eps)[0],
        }
        ibd_alt[p] = parts
        coh = {a: ibd_cohesion(W, l) for a, l in parts.items()}
        picks = {a: select_cluster(l, coh[a], args.min_cluster_size, args.select)
                 for a, l in parts.items()}
        ibd_pick[p] = {a: picks[a][1] for a in parts}
        for a in parts:
            ibd_sel.setdefault(a, {})[p] = [
                members[i] for i in range(len(members)) if picks[a][0][i]]
        lab = parts[args.ibd_cluster]
        keep, chosen, why = picks[args.ibd_cluster]
        sel = [members[i] for i in range(len(members)) if keep[i]]
        selected[p] = sel
        scores[p] = (emb, lab, keep, members)
        sizes, n_out = describe(lab)
        summary.setdefault(p, {})["ibd"] = {
            "n_total": len(members), "algo": args.ibd_cluster,
            "bin_cm": [args.ibd_min_cm, None if not np.isfinite(args.ibd_max_cm)
                       else args.ibd_max_cm],
            "K": len(sizes), "cluster_sizes": sizes, "n_noise": n_out,
            "n_selected": len(sel), "select": args.select,
            "min_cluster_size": args.min_cluster_size,
            "chosen_cluster": chosen, "why": why,
            "cohesion": {a: {str(k): round(v, 3) for k, v in coh[a].items()}
                         for a in parts},
            "alt": {a: describe(l)[0] for a, l in parts.items()},
            "alt_selected": {a: int(picks[a][0].sum()) for a in parts}}
        print(f"[ibd] {p}: n={len(members)} algo={args.ibd_cluster} "
              f"K={len(sizes)} noise={n_out} -> keep {len(sel)} ({why})")
        for a, l in parts.items():
            print(f"        {a:9s} {compo(l, coh[a], picks[a][1])}")
    methods["ibd"] = (selected, scores)

    # ---- SNP: standard popgen recipe, standardize -> PCA -> K-means ----
    selected, scores = {}, {}
    for p in pops:
        members = [s for s in pop_members[p] if s in gidx]
        feat, ev = pca(D[[gidx[s] for s in members]], args.n_pc)
        lab, K, sil = kmeans_pcs(feat, args.kmax, args.sil_thresh, args.seed)
        coh = snp_cohesion(feat, lab)
        keep, chosen, why = select_cluster(lab, coh, args.min_cluster_size, args.select)
        sel = [members[i] for i in range(len(members)) if keep[i]]
        selected[p] = sel
        scores[p] = (feat, lab, keep, members)
        sizes, n_out = describe(lab)
        var_pct = [round(float(v) * 100, 2) for v in ev[:max(6, args.n_pc + 3)]]
        summary.setdefault(p, {})["snp"] = {
            "n_total": len(members), "n_pc": args.n_pc, "K": int(K),
            "algo": "kmeans", "silhouette": round(sil, 3),
            "cluster_sizes": sizes, "n_selected": len(sel), "var_pct": var_pct,
            "select": args.select, "min_cluster_size": args.min_cluster_size,
            "chosen_cluster": chosen, "why": why,
            "cohesion": {str(k): round(v, 3) for k, v in coh.items()}}
        print(f"[snp] {p}: n={len(members)} kmeans K={K} sil={sil:.2f} "
              f"-> keep {len(sel)} ({why})")
        print(f"        {compo(lab, coh, chosen)}")
        print(f"        var%={var_pct}")
    methods["snp"] = (selected, scores)

    # ---- how much does the IBD algorithm choice matter? ----
    print("\n[compare] IBD partitions: adjusted Rand index (full partition) / "
          "Jaccard of the KEPT set")
    algos = ["louvain", "spectral", "dbscan"]
    for p in pops:
        row = []
        for i in range(len(algos)):
            for j in range(i + 1, len(algos)):
                a, b = ibd_alt[p][algos[i]], ibd_alt[p][algos[j]]
                ari = adjusted_rand_score(a, b)
                ka, kb = largest_cluster_mask(a), largest_cluster_mask(b)
                un = (ka | kb).sum()
                jac = (ka & kb).sum() / un if un else 0.0
                row.append(f"{algos[i][:4]}/{algos[j][:4]} {ari:5.2f}/{jac:4.2f}")
        print(f"  {p}: " + "   ".join(row))
        summary[p]["ibd"]["compare"] = row

    # ---- score every partition against the reported 1000G subpopulations ----
    if truth:
        print("\n[truth] agreement with reported 1000G subpopulations (ARI / AMI)")
        print(f"  {'pop':4s} {'#sub':>4s} | " +
              " ".join(f"{a:>13s}" for a in algos + ["snp-kmeans"]))
        for p in pops:
            mem = pop_members[p]
            y = np.array([truth.get(s.split("_")[0], "?") for s in mem])
            cells = [f"{adjusted_rand_score(y, ibd_alt[p][a]):.2f}/"
                     f"{adjusted_mutual_info_score(y, ibd_alt[p][a]):.2f}"
                     for a in algos]
            smem = [s for s in mem if s in gidx]
            ys = np.array([truth.get(s.split("_")[0], "?") for s in smem])
            sl = methods["snp"][1][p][1]
            cells.append(f"{adjusted_rand_score(ys, sl):.2f}/"
                         f"{adjusted_mutual_info_score(ys, sl):.2f}")
            print(f"  {p:4s} {len(set(y)):>4d} | " +
                  " ".join(f"{c:>13s}" for c in cells))
            summary[p]["truth_ari_ami"] = dict(zip(algos + ["snp-kmeans"], cells))

    # ---- write sliced label files ----
    # <out>_ibd.txt / _snp.txt are the headline pair (IBD via --ibd-cluster),
    # plus one file per IBD algorithm so all three can be fitted and compared
    # rather than trusting the choice of algorithm.
    writes = {m: sel for m, (sel, _) in methods.items()}
    writes.update({f"ibd_{a}": sl for a, sl in ibd_sel.items()})
    for method, sel in writes.items():
        out = f"{args.out}_{method}.txt"
        with open(out, "w") as f:
            for p in pops:
                for s in sel[p]:
                    f.write(f"{s} {p}\n")
        sizes = "/".join(str(len(sel[p])) for p in pops)
        print(f"[save] {out}  ({sum(len(sel[p]) for p in pops)} individuals; "
              f"{'/'.join(pops)} = {sizes})")

    # ---- diagnostic PCA scatter ----
    # IBD rows are reordered co-ancestry heatmaps (the matrix has no coordinate
    # space to scatter); the SNP row stays a PC scatter, where genotypes really
    # are a data matrix and the PC axes are meaningful.
    nr = len(algos) + 1
    fig, axes = plt.subplots(nr, len(pops), figsize=(3.4 * len(pops), 3.4 * nr),
                             squeeze=False)
    for r, algo in enumerate(algos):
        for c, p in enumerate(pops):
            W = coanc[p].copy(); np.fill_diagonal(W, 0.0)
            lab = ibd_alt[p][algo]
            k = len(set(lab[lab >= 0]))
            pick = ibd_pick[p][algo]
            n_sel = int((lab == pick).sum()) if pick is not None else 0
            plot_coancestry_blocks(
                axes[r][c], W, lab,
                f"IBD {p} -- {algo} (K={k}, n={len(lab)}, keep {n_sel})",
                chosen=pick)
    for c, p in enumerate(pops):
        ax = axes[nr - 1][c]
        feat, lab, keep, _ = methods["snp"][1][p]
        ax.scatter(feat[~keep, 0], feat[~keep, 1], s=8, c="#bbbbbb", label="dropped")
        ax.scatter(feat[keep, 0], feat[keep, 1], s=8, c="#c0392b", label="kept")
        # Clip the view to the bulk: a single extreme sample (AFR has one at
        # PC1 ~ -140) otherwise squashes every real cluster into a sliver.
        # Percentiles, not MAD: these PC distributions are deliberately
        # multimodal, which inflates the MAD and would zoom OUT past the data.
        off = np.zeros(len(feat), bool)
        for axis, setlim in ((0, ax.set_xlim), (1, ax.set_ylim)):
            lo, hi = np.percentile(feat[:, axis], [1, 99])
            pad = 0.08 * (hi - lo) + 1e-9
            lo, hi = lo - pad, hi + pad
            off |= (feat[:, axis] < lo) | (feat[:, axis] > hi)
            setlim(lo, hi)
        n_off = int(off.sum())
        ax.set_title(f"SNP {p} -- kmeans K={len(set(lab[lab >= 0]))} "
                     f"(keep {keep.sum()}/{len(keep)}"
                     + (f", {n_off} off-view)" if n_off else ")"), fontsize=9)
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    axes[nr - 1][0].legend(fontsize=8, loc="best")
    fig.suptitle(f"IBD: co-ancestry enrichment log2(W / k_i k_j 2m) reordered by "
                 f"cluster, {args.ibd_min_cm:g}+ cM  (red diagonal blocks = real "
                 f"communities)   |   SNP: PC scatter\n"
                 f"green box / red points = cluster kept by --select {args.select} "
                 f"(min size {args.min_cluster_size})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{args.out}_clusters.png", dpi=150, bbox_inches="tight")
    print(f"[save] {args.out}_clusters.png")

    json.dump(summary, open(f"{args.out}_summary.json", "w"), indent=2)
    print(f"[save] {args.out}_summary.json")


if __name__ == "__main__":
    main()
