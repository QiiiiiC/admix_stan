"""Dense-data ADVI refit on the TWO-ADMIXTURE 4-pop topology, HapNe-style bins
(2-20.5 cM uniform 0.5 cM).  Fits every model variant -- Nfixed / Nvarying /
Nepoch / Nsmooth (random-walk smoothness prior) -- and writes, into --out-dir:
  dense_{tag}_spectrum{sfx}.png   IBD spectrum (value)
  dense_{tag}_density{sfx}.png    IBD density (value / cM)
  dense_{tag}_elbo{sfx}.png       ADVI ELBO vs iteration
  dense_{tag}_params{sfx}.json    posterior-mean params + ELBO trajectory

Incremental cache: an existing params json is loaded and only the MISSING models
are (re)fitted, so adding the Nsmooth models on top of a cached 13-model run only
costs the 4 new fits.

    PY=/opt/miniconda3/envs/genetics_env/bin/python
    $PY fit_dense.py --out-dir results/raw
    $PY fit_dense.py --out-dir results/ibd_cluster \
        --allchr-prefix real_4pop_allchr_ibdcluster --chr1-prefix real_chr1_4pop_ibdcluster
"""
import os, sys, json, glob, re, tempfile, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RD = os.path.dirname(os.path.abspath(__file__))          # new_pipeline/real_data
MD = os.path.join(os.path.dirname(RD), "models")         # new_pipeline/models
ADVI_DIR = os.path.join(tempfile.gettempdir(), "fit_dense_advi")
os.makedirs(ADVI_DIR, exist_ok=True)
sys.path.insert(0, RD)
from cmdstanpy import CmdStanModel
import infer_topology as IT
import topologies as T

POPS = ["AFR", "EAS", "EUR", "SAS"]

SPECS = [
    ("IBD-only",         "ibd_model_Nfixed.stan",           "ibd",   False, "fixed"),
    ("IBD-masked",       "ibd_model_Nfixed_masked.stan",    "ibd",   False, "fixed"),
    ("Mixed",            "mixed_model_Nfixed.stan",         "mixed", True,  "fixed"),
    ("Mixed-masked",     "mixed_model_Nfixed_masked.stan",  "mixed", True,  "fixed"),
    ("SNP-only",         "snp_model_Nfixed.stan",           "snp",   False, "fixed"),
    ("IBD-only-Nv",      "ibd_model_Nvarying.stan",         "ibd",   False, "nvary"),
    ("IBD-masked-Nv",    "ibd_model_Nvarying_masked.stan",  "ibd",   False, "nvary"),
    ("Mixed-Nv",         "mixed_model_Nvarying.stan",       "mixed", True,  "nvary"),
    ("Mixed-masked-Nv",  "mixed_model_Nvarying_masked.stan","mixed", True,  "nvary"),
    ("IBD-only-Nep",     "ibd_model_Nepoch.stan",           "ibd",   False, "epoch"),
    ("IBD-masked-Nep",   "ibd_model_Nepoch_masked.stan",    "ibd",   False, "epoch"),
    ("Mixed-Nep",        "mixed_model_Nepoch.stan",         "mixed", True,  "epoch"),
    ("Mixed-masked-Nep", "mixed_model_Nepoch_masked.stan",  "mixed", True,  "epoch"),
    ("IBD-only-Nsm",     "ibd_model_Nsmooth.stan",          "ibd",   False, "smooth"),
    ("IBD-masked-Nsm",   "ibd_model_Nsmooth_masked.stan",   "ibd",   False, "smooth"),
    ("Mixed-Nsm",        "mixed_model_Nsmooth.stan",        "mixed", True,  "smooth"),
    ("Mixed-masked-Nsm", "mixed_model_Nsmooth_masked.stan", "mixed", True,  "smooth"),
]
ORDER = [s[0] for s in SPECS]


def _find_model(sf):
    for d in (MD, RD):                 # Nsmooth live in real_data/, the rest in models/
        p = os.path.join(d, sf)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(sf)


MODELS = {sf: CmdStanModel(stan_file=_find_model(sf)) for _, sf, *_ in SPECS}
FWD = {"fixed": CmdStanModel(stan_file=_find_model("ibd_model_Nfixed.stan")),
       "nvary": CmdStanModel(stan_file=_find_model("ibd_model_Nvarying.stan")),
       "epoch": CmdStanModel(stan_file=_find_model("ibd_model_Nepoch.stan"))}


def smooth_tree_data(dem):
    """Tree structure for the Nsmooth log-Ne random walk, in dem.nodes order
    (1-based node indices).  Each non-root node points at the branch it flows
    into going backwards in time: a single parent (MERGE child) or, for an
    admixture source node, its admixture index (two parents via admixture_map)."""
    names = list(dem.nodes.keys())
    idx = {nm: i + 1 for i, nm in enumerate(names)}
    N = len(names)
    ne_parent = [0] * N
    ne_admix_idx = [0] * N
    ne_start_event = [0] * N
    a_i = 0
    for k, ev in enumerate(dem.ordered_events):
        ev_no = k + 1                                     # matches cumulative_times[ev_no]
        if ev["type"] == "MERGE":
            P = ev["parent"]; c1, c2 = ev["children"]
            ne_start_event[idx[P] - 1] = ev_no
            ne_parent[idx[c1] - 1] = idx[P]
            ne_parent[idx[c2] - 1] = idx[P]
        else:                                             # ADMIXTURE
            a_i += 1
            X = ev["child"]; s1, s2 = ev["parents"]
            ne_start_event[idx[s1] - 1] = ev_no
            ne_start_event[idx[s2] - 1] = ev_no
            ne_admix_idx[idx[X] - 1] = a_i
    return {"ne_parent": ne_parent, "ne_admix_idx": ne_admix_idx,
            "ne_start_event": ne_start_event}


def make_init(kind, kappa, NE_EV, NE_ADM, N_NODES, N_EPOCH):
    init = {"times": [100.0]*NE_EV, "admixture_fractions": [0.5]*NE_ADM}
    if kind == "fixed":
        init["effective_N"] = 15000.0
    elif kind in ("nvary", "smooth"):
        init.update(mu_log=float(np.log(15000)), sigma_log=0.3, Ne_raw=[0.0]*N_NODES)
        if kind == "smooth":
            init["tau"] = 0.3
    elif kind == "epoch":
        init["Ne_epoch"] = [15000.0]*N_EPOCH
    # `kappa` is kept in the signature (callers pass it per-spec) but is now a
    # no-op: the mixed models take w_se at face value and add the two likelihoods
    # directly.  See the SNP-likelihood comment in mixed_model_*.stan.
    return init


def summarize(sv, kind, dem, NE_ADM):
    p = {"kind": kind,
         "times": np.asarray(sv["times"]).mean(0).tolist(),
         "admix": (np.asarray(sv["admixture_fractions"]).reshape(-1, NE_ADM).mean(0).tolist()
                   if NE_ADM else [])}
    if kind == "fixed":
        p["effective_N"] = float(np.mean(sv["effective_N"]))
        p["ne_recent"] = p["ne_deep"] = p["effective_N"]
    elif kind in ("nvary", "smooth"):
        p["Ne_vec"] = np.asarray(sv["Ne"]).mean(0).tolist()
        est = IT.extract_estimates(dem, sv, True)["effective_sizes"]
        p["ne_recent"] = float(np.mean([est[q]["Ne_haploid_mean"] for q in POPS]))
        p["ne_deep"] = est["root"]["Ne_haploid_mean"]
        if "tau" in sv:
            p["tau"] = float(np.mean(sv["tau"]))
    elif kind == "epoch":
        p["Ne_epoch"] = np.asarray(sv["Ne_epoch"]).mean(0).tolist()
        p["ne_recent"], p["ne_deep"] = p["Ne_epoch"][0], p["Ne_epoch"][-1]
    if "kappa_snp" in sv:
        p["kappa"] = float(np.mean(sv["kappa_snp"]))
    return p


def parse_elbo(outdir):
    files = sorted(glob.glob(os.path.join(outdir, "*-stdout.txt")),
                   key=os.path.getmtime)
    if not files:
        return [], []
    it, el = [], []
    for line in open(files[-1], errors="ignore"):
        mobj = re.match(r"\s*(\d+)\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s", line)
        if mobj:
            it.append(int(mobj.group(1))); el.append(float(mobj.group(2)))
    return it, el


def fit_advi(sf, data, init, kind, dem, NE_ADM, tag, name):
    outdir = os.path.join(ADVI_DIR, f"{tag}_{name}")
    os.makedirs(outdir, exist_ok=True)
    for f in glob.glob(os.path.join(outdir, "*")):
        os.remove(f)
    for seed in (1, 7, 13):
        try:
            fit = MODELS[sf].variational(
                data=data, inits=init, seed=seed, iter=5000, eval_elbo=50,
                grad_samples=1, elbo_samples=100, tol_rel_obj=1e-3,
                require_converged=False, output_dir=outdir, draws=200,
                show_console=False)
            sv = fit.stan_variables(mean=False)
            p = summarize(sv, kind, dem, NE_ADM)
            p["elbo_iters"], p["elbo_vals"] = parse_elbo(outdir)
            return p
        except Exception as ex:
            print(f"      [ADVI seed {seed}] {str(ex).splitlines()[-1][:60]}", flush=True)
    return None


def predict(p, ibd_data):
    base = {"times": p["times"], "admixture_fractions": p["admix"]}
    if p["kind"] == "fixed":
        init = {**base, "effective_N": p["effective_N"]}
    elif p["kind"] == "epoch":
        init = {**base, "Ne_epoch": p["Ne_epoch"]}
    else:                                   # nvary or smooth -> reproduce Ne vector
        Ne = np.asarray(p["Ne_vec"]); init = {**base, "mu_log": 0.0, "sigma_log": 1.0, "Ne_raw": np.log(Ne).tolist()}
    fwd_kind = "nvary" if p["kind"] == "smooth" else p["kind"]
    f = FWD[fwd_kind].sample(data=ibd_data, fixed_param=True, iter_sampling=1, iter_warmup=0,
                             adapt_engaged=False, chains=1, inits=init, show_console=False)
    return np.asarray(f.stan_variable("ibd_fraction"))[0]


COL = {"IBD-only":"#08519c","IBD-masked":"#3182bd","Mixed":"#6baed6","Mixed-masked":"#bdd7e7",
       "SNP-only":"#2ca02c",
       "IBD-only-Nv":"#54278f","IBD-masked-Nv":"#807dba","Mixed-Nv":"#9e9ac8","Mixed-masked-Nv":"#cbc9e2",
       "IBD-only-Nep":"#a63603","IBD-masked-Nep":"#e6550d","Mixed-Nep":"#fd8d3c","Mixed-masked-Nep":"#fdd0a2",
       "IBD-only-Nsm":"#01665e","IBD-masked-Nsm":"#35978f","Mixed-Nsm":"#80cdc1","Mixed-masked-Nsm":"#c7eae5"}
LW = {n:(5 if n.startswith("IBD-only") else 3.3 if n.startswith("IBD-masked") else 2 if n.startswith("Mixed") and "masked" not in n else 1.2) for n in ORDER}
LW["SNP-only"]=2.2
LS = {n:("--" if "masked" in n else "-") for n in ORDER}


def spectrum_plot(preds, obs, se, xmid, xlo, names, title, out, density=False, width=None):
    fig, axes = plt.subplots(4, 4, figsize=(19, 14), sharex=True)
    for i in range(4):
        for j in range(4):
            ax = axes[i][j]
            if j < i:
                ax.axis("off"); continue
            for n in names:
                y = preds[n][:, i, j]/width if density else preds[n][:, i, j]
                ax.plot(xmid, y, LS[n], color=COL[n], lw=LW[n], alpha=0.9, label=n)
            o = obs[:, i, j]; s = se[:, i, j]; nz = o > 0
            yo = (o/width) if density else o; so = (s/width) if density else s
            ax.errorbar(xmid[nz], yo[nz], yerr=so[nz], fmt="o", color="k", ms=4, capsize=2, lw=1, label="obs", zorder=10)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"{POPS[i]}-{POPS[j]}", fontsize=10)
            ax.set_xlim(xlo*0.95, 22)
            ax.set_ylim(1e-10,1e-3) if density else ax.set_ylim(1e-9,1e-2)
            if i==3 or j==i: ax.set_xlabel("segment length (cM)")
            if j==i: ax.set_ylabel("IBD fraction/cM" if density else "IBD fraction")
    h,l = axes[0][0].get_legend_handles_labels()
    axes[2][0].legend(h,l,loc="center",fontsize=8,title="blue=Nfixed purple=Nvarying\norange=Nepoch teal=Nsmooth")
    fig.suptitle(title, fontsize=14); fig.tight_layout(rect=[0,0,1,0.97])
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig); print(f"[save] {out}")


def run(tag, prefix, out_dir, sfx="", topology="two_admix", specs=None):
    print("\n"+"#"*70+f"\n### {tag}{sfx}  ->  {out_dir}\n"+"#"*70, flush=True)
    os.makedirs(out_dir, exist_ok=True)
    dem = {"two_admix": T.build_two_admix_4pop,
           "no_admix":  T.build_tree_no_admix_4pop}[topology](POPS)
    specs = specs or SPECS
    order = [s[0] for s in specs]
    ibd_data   = IT.assemble_stan_data(dem, prefix+".npz", prefix+".json", model="ibd",   T_max=IT.DEFAULT_T_MAX)
    mixed_data = IT.assemble_stan_data(dem, prefix+".npz", prefix+".json", model="mixed", T_max=IT.DEFAULT_T_MAX)
    snp_data   = IT.assemble_stan_data(dem, prefix+".npz", prefix+".json", model="snp",   T_max=IT.DEFAULT_T_MAX)
    tree = smooth_tree_data(dem)                 # Nsmooth needs the parent/edge structure
    ibd_data.update(tree); mixed_data.update(tree)
    NE_EV, NE_ADM, N_NODES = ibd_data["n_events"], ibd_data["n_admixture"], ibd_data["n_nodes"]
    N_EPOCH = NE_EV + 1
    DATA = {"ibd": ibd_data, "mixed": mixed_data, "snp": snp_data}
    CACHE = os.path.join(out_dir, f"dense_{tag}_params{sfx}.json")

    fits = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    changed = False
    for name, sf, dk, kap, kind in specs:
        if fits.get(name, {}).get("advi") is not None:
            continue                                   # already cached
        print(f"[{tag}{sfx}] fit {name}", flush=True)
        init = make_init(kind, kap, NE_EV, NE_ADM, N_NODES, N_EPOCH)
        ad = fit_advi(sf, DATA[dk], init, kind, dem, NE_ADM, tag, name)
        fits[name] = {"advi": ad}; changed = True
    if changed:
        json.dump(fits, open(CACHE, "w"), indent=2); print(f"[save] {CACHE}")

    npz = np.load(prefix+".npz"); obs = npz["ibd_mean"]; se = np.sqrt(npz["ibd_var"])
    binL = np.array(ibd_data["bin_length"]); xmid = np.sqrt(binL[:,0]*binL[:,1]); xlo = float(binL[:,0].min())
    width = binL[:,1]-binL[:,0]

    got = {n: fits[n]["advi"] for n in order if fits.get(n, {}).get("advi")}
    names = [n for n in order if n in got]
    preds = {n: predict(got[n], ibd_data) for n in names}
    spectrum_plot(preds, obs, se, xmid, xlo, names,
                  f"Dense {tag}{sfx} IBD spectrum (ADVI, {topology}, HapNe bins)",
                  os.path.join(out_dir, f"dense_{tag}_spectrum{sfx}.png"))
    spectrum_plot(preds, obs, se, xmid, xlo, names,
                  f"Dense {tag}{sfx} IBD density (ADVI, {topology}, HapNe bins)",
                  os.path.join(out_dir, f"dense_{tag}_density{sfx}.png"), density=True, width=width)

    # ---- ELBO vs iteration (ADVI), small multiples ----
    ncol = min(5, len(order))
    nrow = (len(order) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 3.2*nrow), squeeze=False)
    axes = axes.ravel()
    for k, n in enumerate(order):
        ax = axes[k]; ad = fits.get(n, {}).get("advi")
        if ad and ad.get("elbo_iters"):
            ax.plot(ad["elbo_iters"], ad["elbo_vals"], "-", color=COL[n], lw=1.5)
            ax.set_title(n, fontsize=9); ax.set_xlabel("iter"); ax.set_ylabel("ELBO")
            ax.set_yscale("symlog")
        else:
            ax.text(0.5,0.5,f"{n}\n(no ADVI)",ha="center",va="center"); ax.axis("off")
    for k in range(len(order), len(axes)): axes[k].axis("off")
    fig.suptitle(f"Dense {tag}{sfx}: ADVI ELBO vs iteration ({topology})", fontsize=14); fig.tight_layout(rect=[0,0,1,0.97])
    fig.savefig(os.path.join(out_dir, f"dense_{tag}_elbo{sfx}.png"), dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[save] {os.path.join(out_dir, f'dense_{tag}_elbo{sfx}.png')}")

    # ---- ADVI summary table ----
    print(f"\n=== {tag}{sfx}: ADVI ({topology}, HapNe bins) ===")
    print(f"{'model':<18}{'Ne_recent':>12}{'Ne_deep':>12}{'admix0':>10}{'admix1':>10}{'tau':>7}{'kappa':>10}")
    for n in order:
        ad = fits.get(n, {}).get("advi")
        def g(d,k,default=float('nan')):
            return d.get(k,default) if d else float('nan')
        a0=(ad['admix'][0] if ad and ad.get('admix') else float('nan'))
        a1=(ad['admix'][1] if ad and ad.get('admix') and len(ad['admix'])>1 else float('nan'))
        print(f"{n:<18}{g(ad,'ne_recent'):>12.0f}{g(ad,'ne_deep'):>12.0f}{a0:>10.3f}{a1:>10.3f}{g(ad,'tau'):>7.2f}{g(ad,'kappa'):>10.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="ADVI two-admix refit incl. Nsmooth models")
    ap.add_argument("--out-dir", default=RD, help="Directory for dense_* outputs.")
    ap.add_argument("--suffix", default="", help="Optional filename suffix.")
    ap.add_argument("--allchr-prefix", default=os.path.join(RD, "real_4pop_allchr"))
    ap.add_argument("--chr1-prefix", default=os.path.join(RD, "real_chr1_4pop"))
    ap.add_argument("--topology", choices=["two_admix", "no_admix"],
                    default="two_admix",
                    help="no_admix = the plain ((EUR,SAS),EAS),AFR) tree.")
    ap.add_argument("--models", nargs="+", default=None, metavar="NAME",
                    help=f"Subset of models to fit. Default all. Choices: {ORDER}")
    ap.add_argument("--tags", nargs="+", default=["allchr", "chr1"],
                    help="Which datasets to run (allchr and/or chr1).")
    args = ap.parse_args()
    sfx = f"_{args.suffix}" if args.suffix else ""
    specs = SPECS
    if args.models:
        unknown = [m for m in args.models if m not in ORDER]
        if unknown:
            ap.error(f"unknown model(s) {unknown}; choose from {ORDER}")
        specs = [s for s in SPECS if s[0] in set(args.models)]
    prefixes = {"allchr": args.allchr_prefix, "chr1": args.chr1_prefix}
    for tag in args.tags:
        run(tag, prefixes[tag], args.out_dir, sfx, args.topology, specs)
    print("\nALL DONE")
