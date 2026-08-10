"""Screen candidate leaves for admixture with the f3 statistic, before fitting anything.

A tree leaf must have non-negative branch drift.  f3(X; A, B) IS that branch drift
when the tree is ((X,A),B), so a significantly NEGATIVE f3(X; A, B) says X cannot
be a leaf of the tree relating X, A and B -- it is a mixture of lineages related
to A and to B (Patterson et al. 2012, Genetics 192:1065).  This is invisible on
two leaves: with L populations a double-centered covariance has L(L-1)/2 free
numbers and an unrooted tree has 2L-3 branches, so L=2 gives 1 and 1 -- you
recover only x_A + x_B and can never test either separately.  L=3 gives 3 and 3.

Estimators are the unbiased ones, which subtract each population's OWN per-SNP
heterozygosity:

    f2(A,B) = mean_a [ (pA - pB)^2 - pA(1-pA)/(nA-1) - pB(1-pB)/(nB-1) ]
    f3(C;A,B) = mean_a [ (pC - pA)(pC - pB) - pC(1-pC)/(nC-1) ]

NOT the TreeMix-normalised covariance: its bias correction is a uniform scalar
1/n_haploid, which contributes exactly -1/n_target to f3 and so can manufacture
the sign it is being used to test.  Errors are a delete-one-chromosome jackknife
(21 blocks, independent by construction).
"""
import os, sys, glob, json, argparse
import numpy as np

RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RD)
import generate_stan_data as G

ap = argparse.ArgumentParser()
ap.add_argument("--pops", nargs="+", required=True, help="candidate leaves")
ap.add_argument("--outgroup", required=True)
ap.add_argument("--labels", default="high_cov/samples/labels_subpop.txt")
ap.add_argument("--vcf-glob", default="merged_pruned/vcf/merged_all.chr[!9]*_pruned.vcf.gz")
ap.add_argument("--min-maf", type=float, default=0.05)
ap.add_argument("--cache", default="clustering/results_screen/f3_screen_freq.npz")
a = ap.parse_args()
os.chdir(RD)

POPS = list(a.pops) + [a.outgroup]
p2i = {p: i for i, p in enumerate(POPS)}
s2p = {}
for line in open(a.labels):
    if line.strip():
        s, p = line.split()[:2]
        if p in p2i:
            s2p[s] = p2i[p]
n_ind = np.bincount([s2p[s] for s in s2p], minlength=len(POPS))
n_hap = (G.PLOIDY * n_ind).astype(float)
print(f"[labels] " + "  ".join(f"{p}:{n}" for p, n in zip(POPS, n_ind)))

if os.path.exists(a.cache):
    z = np.load(a.cache, allow_pickle=True)
    F, chrom, cpops = z["F"].astype(np.float64), z["chrom"], list(z["pops"])
    assert cpops == POPS, f"cache is for {cpops}, not {POPS}; delete it"
    print(f"[cache] {a.cache}: {len(F):,} SNPs")
else:
    rows, chs = [], []
    for v in sorted(glob.glob(a.vcf_glob)):
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
        print(f"[vcf] {os.path.basename(v)} -> {len(rows):,} sites", flush=True)
    F = np.asarray(rows); chrom = np.asarray(chs)
    mu = np.nanmean(F, axis=1)
    ok = (mu > a.min_maf) & (mu < 1 - a.min_maf) & ~np.isnan(F).any(1)
    F, chrom = F[ok], chrom[ok]
    print(f"[snp] {ok.sum():,}/{len(mu):,} SNPs kept (min_maf={a.min_maf})")
    np.savez_compressed(a.cache, F=F.astype(np.float32), chrom=chrom, pops=np.array(POPS))

h = F * (1 - F) / (n_hap - 1)[None, :]
grp = [np.where(chrom == c)[0] for c in sorted(set(chrom))]
alln = np.arange(len(F))


def jk(fn):
    full = fn(alln)
    lo = np.array([fn(np.setdiff1d(alln, g)) for g in grp])
    m = len(grp)
    return full, np.sqrt((m - 1) / m * ((lo - lo.mean()) ** 2).sum())


f3 = lambda I, c, x, y: np.mean((F[I, c] - F[I, x]) * (F[I, c] - F[I, y]) - h[I, c])
f2 = lambda I, i, j: np.mean((F[I, i] - F[I, j]) ** 2 - h[I, i] - h[I, j])

O = len(POPS) - 1
C = list(range(O))
print(f"\nf3(X; Y, {a.outgroup}) -- X's own branch drift under the tree ((X,Y),{a.outgroup}).")
print(f"NEGATIVE and |Z| > 3  =>  X is admixed and cannot be a clean leaf.\n")
print(f"{'X':>6} {'Y':>6} {'f3':>11} {'SE':>10} {'Z':>8}   verdict")
res = {}
worst = {p: (1e9, None) for p in POPS[:O]}
for x in C:
    for y in C:
        if x == y:
            continue
        v, s = jk(lambda I, c=x, i=y: f3(I, c, i, O))
        z = v / s
        res[f"f3({POPS[x]};{POPS[y]},{a.outgroup})"] = {"f3": v, "se": s, "Z": z}
        flag = "ADMIXED" if z < -3 else ("ok" if z > 3 else "ambiguous")
        print(f"{POPS[x]:>6} {POPS[y]:>6} {v:>+11.6f} {s:>10.6f} {z:>+8.2f}   {flag}")
        if z < worst[POPS[x]][0]:
            worst[POPS[x]] = (z, POPS[y])

print(f"\n{'population':>12} {'worst Z over all Y':>20}   status")
clean = []
for p in POPS[:O]:
    z, y = worst[p]
    st = "ADMIXED" if z < -3 else "clean"
    if z >= -3:
        clean.append(p)
    print(f"{p:>12} {z:>+13.2f} (Y={y:<4}) {'':>2} {st}")

print(f"\npairwise f2 among the clean candidates (smaller = more recent split):")
for i in C:
    for j in C:
        if i < j and POPS[i] in clean and POPS[j] in clean:
            v, s = jk(lambda I, a_=i, b_=j: f2(I, a_, b_))
            print(f"   f2({POPS[i]},{POPS[j]}) = {v:.6f}  SE {s:.6f}")

os.makedirs("clustering/results_screen", exist_ok=True)
json.dump(res, open("clustering/results_screen/f3_screen.json", "w"), indent=2)
print(f"\n[save] clustering/results_screen/f3_screen.json")
