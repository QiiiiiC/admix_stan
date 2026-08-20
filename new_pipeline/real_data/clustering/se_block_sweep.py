"""Measure kappa_snp directly, without the demographic model.

kappa_snp exists because the block SE the TreeMix estimator reports may be too
small: `resample_snp_covariance` splits the concatenated SNP list into blocks of
`se_block_size` (default 500) and returns sd(block means)/sqrt(B).  That is only
the true sampling SE if the blocks are independent.  Residual LD makes nearby
blocks correlated, and correlated blocks make sd(block means)/sqrt(B) an
UNDER-estimate -- exactly the deficit kappa is meant to absorb.

But that deficit is a property of the DATA, not of the demography, so it can be
measured on its own: enlarge the blocks until the SE stops growing.  Once a block
is much longer than the LD correlation length the estimate must plateau, and

    kappa = SE(plateau) / SE(500 SNPs)

is then a plug-in constant, learned from data and held fixed in the fit -- no
gradient, no prior, no degeneracy.  The chromosome-level split is the reference
point: LD does not cross chromosomes, so those 22 blocks are independent by
construction and their SE is the honest one.

Writes <prefix>_freq.npz (dev, het, chrom) so the sweep and any later refit reuse
one pass over the VCFs.
"""
import os, sys, glob, json, argparse
import numpy as np

RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RD)
import generate_stan_data as G

ap = argparse.ArgumentParser()
ap.add_argument("--labels", required=True)
ap.add_argument("--pop-order", nargs="+", required=True)
ap.add_argument("--vcf-glob", default="merged_pruned/vcf/merged_all.chr[!9]*_pruned.vcf.gz")
ap.add_argument("--min-maf", type=float, default=0.05)
ap.add_argument("--out", required=True)
a = ap.parse_args()

os.chdir(RD)
POPS = a.pop_order
p2i = {p: i for i, p in enumerate(POPS)}
s2p = {}
for line in open(a.labels):
    if line.strip():
        s, p = line.split()[:2]
        if p in p2i:
            s2p[s] = p2i[p]
n_ind = np.bincount([s2p[s] for s in s2p], minlength=len(POPS))
n_haps = (G.PLOIDY * n_ind).astype(float)
print(f"[labels] {dict(zip(POPS, n_ind))}  haploid n = {n_haps}")

cache = a.out + "_freq.npz"
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True)
    dev, het, chrom = z["dev"], z["het"], z["chrom"]
    print(f"[cache] {cache}: {len(het):,} SNPs")
else:
    vcfs = sorted(glob.glob(a.vcf_glob))
    print(f"[vcf] {len(vcfs)} files")
    rows, chs = [], []
    for v in vcfs:
        names = None
        with G.open_any(v) as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    names = line.rstrip("\n").split("\t")[9:]
                    cols = [[] for _ in POPS]
                    for c, s in enumerate(names):
                        if G.base_id(s) in s2p:
                            cols[s2p[G.base_id(s)]].append(c)
                    cols = [np.asarray(c, int) for c in cols]
                    na = np.array([G.PLOIDY * len(c) for c in cols], float)
                    continue
                parts = line.rstrip("\n").split("\t")
                gts = np.asarray(parts[9:], dtype="U3")
                alt = ((gts == "1|0").astype(np.int8) + (gts == "0|1").astype(np.int8)
                       + 2 * (gts == "1|1").astype(np.int8))
                rows.append([alt[c].sum() / n for c, n in zip(cols, na)])
                chs.append(parts[0])
        print(f"[vcf] {os.path.basename(v)} -> {len(rows):,} cumulative sites", flush=True)
    F = np.asarray(rows)
    chrom = np.asarray(chs)
    mu = np.nanmean(F, axis=1)
    ok = (mu > a.min_maf) & (mu < 1 - a.min_maf) & ~np.isnan(mu) & ~np.isnan(F).any(1)
    dev, het, chrom = F[ok] - mu[ok, None], mu[ok] * (1 - mu[ok]), chrom[ok]
    print(f"[snp] {ok.sum():,}/{len(mu):,} SNPs kept (min_maf={a.min_maf})")
    # Keep the raw per-population frequencies too: the unbiased f2 estimators
    # need each population's OWN heterozygosity, which `dev` and the pooled `het`
    # (built from the cross-population mean) have already thrown away.
    np.savez_compressed(cache, dev=dev, het=het, chrom=chrom, F=F[ok].astype(np.float32))
    print(f"[save] {cache}")

P = dev.shape[1]
w_hat = (dev.T @ dev) / het.sum()
bias = np.diag(1.0 / n_haps)
rm = bias.mean(1, keepdims=True)
w_hat = w_hat - (bias - rm - rm.T + bias.mean())


def se_from_groups(groups):
    """SE of the mean over a set of index groups, TreeMix's estimator."""
    wb = np.stack([(dev[g].T @ dev[g]) / het[g].sum() for g in groups])
    B = len(groups)
    return np.sqrt(((wb - wb.mean(0)) ** 2).sum(0) / (B * (B - 1))), B


print(f"\n{'blocks':>22} {'B':>5} " + " ".join(f"{POPS[i]+'-'+POPS[j]:>12}"
      for i in range(P) for j in range(i, P)))
rowfmt = lambda tag, B, se: (f"{tag:>22} {B:>5} " +
    " ".join(f"{se[i,j]:12.3e}" for i in range(P) for j in range(i, P)))

res = {}
J = len(het)
base = None
for bs in (100, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000):
    B = J // bs
    if B < 3:
        continue
    idx = np.arange(B * bs).reshape(B, bs)
    se, B = se_from_groups(list(idx))
    res[f"snp{bs}"] = se.tolist()
    if bs == 500:
        base = se
    print(rowfmt(f"{bs} SNPs", B, se))

# LD cannot cross a chromosome boundary, so these blocks are independent by
# construction -- this is the SE every other row should be compared against.
cg = [np.where(chrom == c)[0] for c in sorted(set(chrom))]
se_chr, Bc = se_from_groups(cg)
res["chrom"] = se_chr.tolist()
print(rowfmt("per chromosome", Bc, se_chr))

# Chromosome ARMS: 39 independent units instead of 22, the HapNe convention.
# More blocks means a better-conditioned SE with the same independence guarantee.
print()
print(f"kappa = SE(per-chromosome) / SE(500-SNP blocks), per cell:")
K = se_chr / base
for i in range(P):
    for j in range(i, P):
        print(f"   {POPS[i]}-{POPS[j]:<5} {K[i,j]:6.2f}")
print(f"\n   overall (rms over the {P*(P+1)//2} cells): {np.sqrt((K[np.triu_indices(P)]**2).mean()):.2f}")

json.dump({"pops": POPS, "w_hat": w_hat.tolist(), "se": res,
           "kappa_chrom_over_500": K.tolist()},
          open(a.out + "_se_sweep.json", "w"), indent=2)
print(f"[save] {a.out}_se_sweep.json")
