"""Three-leaf fit: Pathfinder vs NUTS, with the composite likelihood decomposed.

Three leaves is the smallest graph on which the SNP branch drifts are separately
identifiable (see topologies.build_tree_3pop), and the leaves are screened for
admixture with f3 first (clustering/f3_screen.py), because a negative f3 means the
population cannot be a tree leaf at all -- and that is invisible on two leaves.

kappa_snp is gone.  The IBD and SNP log-densities are simply added, which is the
independent-product composite likelihood; w_se is now computed on 16 Mb blocks
(--block-size-mb), matching TreeMix's own physical-length criterion, so the
correct relative weight is 1 by construction rather than by fiat.

The decomposition comes from `generated quantities` inside the Stan models rather
than a Python re-implementation, so it cannot drift away from what was actually
sampled, and it is evaluated PER DRAW rather than at the posterior mean -- see
`report_decomposition` for why that distinction matters here.
"""
import os, sys, json, glob, time, tempfile, argparse
import numpy as np

RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RD)
from cmdstanpy import CmdStanModel
import infer_topology as IT
import topologies as T
import fit_dense as FD

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="yri_esn_fin")
ap.add_argument("--pops", nargs=3, default=["YRI", "ESN", "FIN"])
ap.add_argument("--engines", nargs="+", default=["pathfinder"],
                choices=["pathfinder", "nuts"],
                help="NUTS is off by default: on the three-leaf graph with kappa\nremoved it reproduced Pathfinder to three digits (t_sis 116 vs 117, every Ne\nwithin 4%) at 2453 s against 13 s. The two-leaf disagreement that motivated\nrunning both was created by kappa_snp, which is gone.")
a = ap.parse_args()

OUT = os.path.join(RD, "clustering", f"results_3pop_{a.tag}")
# Stan inputs live inside the result folder they belong to, so a configuration is
# one self-contained directory: stan_data/ + fits.json + report.md + figures.
PRE = os.path.join(OUT, "stan_data", f"stan_3pop_{a.tag}")
os.makedirs(os.path.dirname(PRE), exist_ok=True)
WORK = os.path.join(tempfile.gettempdir(), f"fit_3pop_{a.tag}")

SPECS = [("SNP-only-Nsm", "snp_model_Nsmooth.stan",   "snp"),
         ("IBD-only-Nsm", "ibd_model_Nsmooth.stan",   "ibd"),
         ("Mixed-Nsm",    "mixed_model_Nsmooth.stan", "mixed")]
ENGINES = a.engines

dem = T.build_tree_3pop(a.pops)
DATA = {k: IT.assemble_stan_data(dem, PRE + ".npz", PRE + ".json", model=k,
                                 T_max=IT.DEFAULT_T_MAX)
        for k in ("ibd", "mixed", "snp")}
DATA_ANY = DATA["mixed"]
for d in DATA.values():
    d.update(FD.smooth_tree_data(dem))
NODES = list(dem.nodes.keys())
N_NODES, NE_EV, NE_ADM = (DATA["ibd"]["n_nodes"], DATA["ibd"]["n_events"],
                          DATA["ibd"]["n_admixture"])
npz = np.load(PRE + ".npz")
WH, WS = npz["w_hat"], npz["w_se"]
print(f"[graph] leaves {a.pops} -> nodes {NODES}, {NE_EV} events")


def unpack(draws, cols):
    out = {}
    for j, c in enumerate(cols):
        if c.endswith("__") and c != "lp__":
            continue
        if "[" in c:
            nm, idx = c[:-1].split("[")
            if "," in idx:
                continue
            out.setdefault(nm, {})[int(idx)] = draws[:, j]
        else:
            out[c] = draws[:, j]
    return {k: (np.column_stack([v[i] for i in sorted(v)]) if isinstance(v, dict) else v)
            for k, v in out.items()}


def run(name, sf, dk, engine):
    model = CmdStanModel(stan_file=FD._find_model(sf))
    init = FD.make_init("smooth", False, NE_EV, NE_ADM, N_NODES, NE_EV + 1)
    od = os.path.join(WORK, f"{name}_{engine}")
    os.makedirs(od, exist_ok=True)
    for f in glob.glob(os.path.join(od, "*")):
        os.remove(f)
    t0 = time.time()
    extra = {}
    if engine == "pathfinder":
        fit = model.pathfinder(data=DATA[dk], inits=init, seed=1, output_dir=od,
                               num_paths=8, draws=4000, num_single_draws=1000,
                               psis_resample=True, show_console=False)
        sv = fit.stan_variables()
    else:
        fit = model.sample(data=DATA[dk], inits=init, seed=1, chains=4,
                           parallel_chains=4, iter_warmup=1500, iter_sampling=1500,
                           adapt_delta=0.95, max_treedepth=12, output_dir=od,
                           show_console=False, show_progress=False)
        sv = fit.stan_variables()
        s = fit.summary()
        extra = {"max_rhat": float(s["R_hat"].max(skipna=True)),
                 "min_ess": float(s["ESS_bulk"].min(skipna=True)),
                 "divergences": fit.diagnose().count("divergent")}
        sv["lp__"] = fit.method_variables()["lp__"].T.ravel()
    sv = {k: np.asarray(v) for k, v in sv.items()}
    r = {"secs": time.time() - t0, "engine": engine, **extra}
    t = np.asarray(sv["times"]).reshape(len(sv["times"]), -1)
    Ne = np.asarray(sv["Ne"])
    r["times"] = t.mean(0).tolist()
    r["Ne"] = Ne.mean(0).tolist()
    r["Ne_sd_log"] = np.log(Ne).std(0).tolist()
    r["tau"] = float(np.mean(sv["tau"]))
    # SNP branch drift of each sister leaf: x_i = (time to the sister node)/Ne_i.
    # Their RATIO estimates Ne_1/Ne_0 with no reference to t -- the t-free
    # comparison that two leaves cannot support.
    t_sis = t[:, 0]
    x = np.column_stack([t_sis / Ne[:, 0], t_sis / Ne[:, 1]])
    r["x_drift"] = x.mean(0).tolist()
    r["Ne_ratio_10"] = float(np.mean(Ne[:, 1] / Ne[:, 0]))
    for k in ("lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp", "n_ibd_obs", "n_snp_obs"):
        if k in sv:
            v = np.asarray(sv[k]).ravel().astype(float)
            r[k] = {"mean": float(v.mean()), "sd": float(v.std()),
                    "q05": float(np.quantile(v, .05)), "q95": float(np.quantile(v, .95))}
    # posterior-mean predicted spectrum / covariance, for the fit plots
    if "ibd_fraction" in sv:
        r["ibd_pred"] = np.asarray(sv["ibd_fraction"]).mean(0).tolist()
    if "W_centered" in sv:
        r["W_pred"] = np.asarray(sv["W_centered"]).mean(0).tolist()
    return r


res = {}
for name, sf, dk in SPECS:
    for eng in ENGINES:
        key = f"{name}|{eng}"
        print(f"\n[fit] {key}", flush=True)
        try:
            res[key] = run(name, sf, dk, eng)
            print(f"      {res[key]['secs']:.0f}s", flush=True)
        except Exception as ex:
            print(f"      FAILED: {str(ex).splitlines()[-1][:100]}", flush=True)
            res[key] = None

json.dump({"pops": a.pops, "nodes": NODES, "w_hat": WH.tolist(), "w_se": WS.tolist(),
           "fits": res}, open(os.path.join(OUT, "fits.json"), "w"), indent=2)

# Observed branch drifts straight from w_hat, for the model to be judged against.
# f2 is invariant to the double-centering, so these are exactly the same numbers
# the SNP likelihood sees -- no re-normalisation needed, unlike the Patterson-scale
# f3 used for the admixture screen.
f2o = lambda i, j: WH[i, i] + WH[j, j] - 2 * WH[i, j]
X_OBS = [(f2o(0, 1) + f2o(0, 2) - f2o(1, 2)) / 2,
         (f2o(0, 1) + f2o(1, 2) - f2o(0, 2)) / 2]

P = a.pops
print("\n" + "=" * 104)
print(f"THREE-LEAF FIT  (({P[0]},{P[1]}),{P[2]})   nodes {NODES}")
print("=" * 104)
print(f"{'model':14s} {'engine':11s} {'t_sis':>7s} {'t_root':>8s} "
      + " ".join(f"{'Ne_'+p:>10s}" for p in P)
      + f" {'x_'+P[0]:>9s} {'x_'+P[1]:>9s} {'ratio':>7s} {'tau':>5s} {'Rhat':>5s} {'secs':>6s}")
for name, *_ in SPECS:
    for eng in ENGINES:
        f = res.get(f"{name}|{eng}")
        if not f:
            print(f"{name:14s} {eng:11s}  (failed)"); continue
        tt = f["times"]
        print(f"{name:14s} {eng:11s} {tt[0]:7.0f} {sum(tt):8.0f} "
              + " ".join(f"{f['Ne'][i]:10,.0f}" for i in range(3))
              + f" {f['x_drift'][0]:9.5f} {f['x_drift'][1]:9.5f} "
              f"{f['Ne_ratio_10']:7.2f} {f['tau']:5.2f} "
              f"{f.get('max_rhat', float('nan')):5.2f} {f['secs']:6.0f}")
print(f"\n{'OBSERVED (from w_hat)':26s} {'':7s} {'':8s} {'':32s} "
      f"{X_OBS[0]:9.5f} {X_OBS[1]:9.5f} {'':7s}"
      f"   <- SNP branch drifts the fit must reproduce")


def report_decomposition():
    """How much of the composite log-likelihood does each data type supply?

    Evaluated per draw, not at the posterior mean.  Plugging the posterior mean
    into a nonlinear log-density gives the value AT a point, which by Jensen is
    an optimistic stand-in for the mean of the log-density; worse, on a ridged
    posterior (and the SNP term makes one -- only t/Ne is identified) the
    componentwise mean can sit off the ridge entirely, so it is not even a
    representative point.  Since the draws are already in hand, the honest
    version costs nothing.

    Raw log-density LEVELS are not comparable between the two components either:
    they have different numbers of terms and different units.  chi2/n is, which
    is why it is reported alongside -- it is the mean squared standardized
    residual, so 1.0 means 'fits within its own stated error bars'.
    """
    print("\n" + "=" * 104)
    print("COMPOSITE LIKELIHOOD DECOMPOSITION (per draw; posterior mean +/- sd)")
    print("=" * 104)
    print(f"{'model':14s} {'engine':11s} {'lp_ibd':>22s} {'lp_snp':>22s} "
          f"{'n_ibd':>6s} {'n_snp':>6s} {'chi2/n ibd':>11s} {'chi2/n snp':>11s}")
    for name, *_ in SPECS:
        for eng in ENGINES:
            f = res.get(f"{name}|{eng}")
            if not f:
                continue
            g = lambda k: (f"{f[k]['mean']:+.1f}+-{f[k]['sd']:.1f}" if k in f else "-")
            n = lambda k: (f"{f[k]['mean']:.0f}" if k in f else "-")
            rr = lambda c, m: (f"{f[c]['mean']/f[m]['mean']:.2f}"
                               if c in f and m in f and f[m]['mean'] else "-")
            print(f"{name:14s} {eng:11s} {g('lp_ibd'):>22s} {g('lp_snp'):>22s} "
                  f"{n('n_ibd_obs'):>6s} {n('n_snp_obs'):>6s} "
                  f"{rr('chi2_ibd','n_ibd_obs'):>11s} {rr('chi2_snp','n_snp_obs'):>11s}")
    print("\nchi2/n is the mean squared standardized residual: 1.0 = fits inside its own\n"
          "stated SEs, >>1 = that data type is being fitted badly, <<1 = its SEs are too big.\n"
          "Compare chi2/n across components; do NOT compare lp levels across components.")


report_decomposition()


# ============================ figures + report table ============================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IBD_HAT = np.asarray(DATA["ibd"]["ibd_hat"])
IBD_SE = np.asarray(DATA["ibd"]["ibd_se"])
BINS = np.asarray(DATA["ibd"].get("bin_centers",
                  np.arange(2.25, 2.25 + 0.5 * IBD_HAT.shape[0], 0.5)))[:IBD_HAT.shape[0]]
PAIRS = [(i, j) for i in range(3) for j in range(i, 3)]
PLBL = [f"{P[i]}-{P[j]}" for i, j in PAIRS]
ENG0 = ENGINES[0]


def fig_spectrum():
    """Observed vs fitted IBD length spectrum, and standardized residuals.

    This is where a constant-per-branch Ne fails visibly: a population whose Ne has
    changed produces a MONOTONE TILT in the residuals across segment length (long
    segments = recent coalescence, short = older), which no single Ne can absorb."""
    fits = [(n, res[f"{n}|{ENG0}"]) for n, *_ in SPECS
            if res.get(f"{n}|{ENG0}") and "ibd_pred" in res[f"{n}|{ENG0}"]]
    if not fits:
        return
    fig, ax = plt.subplots(2, 6, figsize=(23, 7.2), sharex=True)
    for k, (i, j) in enumerate(PAIRS):
        obs, se = IBD_HAT[:, i, j], IBD_SE[:, i, j]
        m = obs > 0
        ax[0, k].errorbar(BINS[m], obs[m], yerr=se[m], fmt="o", ms=3, lw=.8,
                          color="k", label="observed", zorder=3)
        for c, (nm, f) in zip(["#d62728", "#1f77b4", "#2ca02c"], fits):
            pr = np.asarray(f["ibd_pred"])[:, i, j]
            ax[0, k].plot(BINS, pr, lw=1.6, color=c, label=nm)
            r = np.where(m, (obs - pr) / np.where(se > 0, se, 1), np.nan)
            ax[1, k].plot(BINS, r, ".-", ms=3, lw=1, color=c)
        ax[0, k].set_yscale("log"); ax[0, k].set_title(PLBL[k], fontsize=11)
        ax[1, k].axhline(0, color="k", lw=.8)
        for y in (-2, 2):
            ax[1, k].axhline(y, color="grey", lw=.6, ls=":")
        ax[1, k].set_xlabel("segment length (cM)")
        chi = np.nanmean(np.where(m, ((obs - np.asarray(fits[-1][1]["ibd_pred"])[:, i, j])
                                      / np.where(se > 0, se, 1)) ** 2, np.nan))
        ax[1, k].text(.97, .04, f"chi2/n = {chi:.1f}", transform=ax[1, k].transAxes,
                      ha="right", fontsize=9,
                      color=("firebrick" if chi > 8 else "black"))
    ax[0, 0].set_ylabel("IBD fraction"); ax[1, 0].set_ylabel("(obs - fit)/SE")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle(f"IBD length-spectrum fit  (({P[0]},{P[1]}),{P[2]})   {ENG0}", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "spectrum_fit.png"), dpi=130)
    print(f"[plot] {OUT}/spectrum_fit.png")


def fig_likelihood():
    """Who supplies the composite likelihood, and who is fitted well."""
    names = [n for n, *_ in SPECS if res.get(f"{n}|{ENG0}")]
    fits = [res[f"{n}|{ENG0}"] for n in names]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    x = np.arange(len(names)); w = .38
    g = lambda f, k: (f[k]["mean"] if k in f else np.nan)
    ax[0].bar(x - w/2, [g(f, "lp_ibd") for f in fits], w, label="IBD", color="#1f77b4")
    ax[0].bar(x + w/2, [g(f, "lp_snp") for f in fits], w, label="SNP", color="#d62728")
    ax[0].set_ylabel("log-likelihood contribution"); ax[0].legend(fontsize=9)
    ax[0].set_title("lp by component\n(levels NOT comparable across components)", fontsize=10)
    n_ibd = np.nanmax([g(f, "n_ibd_obs") for f in fits])
    n_snp = np.nanmax([g(f, "n_snp_obs") for f in fits])
    ax[1].bar(x - w/2, [g(f, "chi2_ibd") / n_ibd for f in fits], w, color="#1f77b4")
    ax[1].bar(x + w/2, [g(f, "chi2_snp") / n_snp for f in fits], w, color="#d62728")
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].set_yscale("log"); ax[1].set_ylabel("chi2 / n  (1.0 = fits its own SEs)")
    ax[1].set_title(f"fit quality per component\nIBD n={n_ibd:.0f}, SNP n={n_snp:.0f}", fontsize=10)
    ax[2].bar([0, 1], [n_ibd, n_snp], .5, color=["#1f77b4", "#d62728"])
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["IBD", "SNP"])
    ax[2].set_ylabel("number of likelihood terms")
    ax[2].set_title(f"term count: IBD outnumbers SNP {n_ibd/max(n_snp,1):.0f}:1", fontsize=10)
    for k in (0, 1):
        ax[k].set_xticks(x); ax[k].set_xticklabels(names, rotation=20, ha="right", fontsize=9)
    fig.suptitle(f"Composite-likelihood decomposition  (({P[0]},{P[1]}),{P[2]})", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "likelihood_split.png"), dpi=130)
    print(f"[plot] {OUT}/likelihood_split.png")


def fig_drift():
    """Fitted vs observed SNP branch drift -- the t-free comparison the outgroup buys."""
    names = [n for n, *_ in SPECS if res.get(f"{n}|{ENG0}")]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(2); w = .8 / (len(names) + 1)
    ax[0].bar(x - .4 + w/2, X_OBS, w, label="observed (w_hat)", color="k")
    for q, nm in enumerate(names):
        ax[0].bar(x - .4 + w*(q + 1.5), res[f"{nm}|{ENG0}"]["x_drift"], w, label=nm)
    ax[0].axhline(0, color="grey", lw=.8)
    ax[0].set_xticks(x); ax[0].set_xticklabels([f"x_{P[0]}", f"x_{P[1]}"])
    ax[0].set_ylabel("branch drift  t / Ne"); ax[0].legend(fontsize=8)
    ax[0].set_title("SNP branch drift: fitted vs observed", fontsize=10)
    for q, nm in enumerate(names):
        f = res[f"{nm}|{ENG0}"]
        ax[1].bar(q, f["Ne"][0], .35, color="#1f77b4")
        ax[1].bar(q + .35, f["Ne"][1], .35, color="#ff7f0e")
        ax[1].bar(q + .70, f["Ne"][2], .35, color="#2ca02c")
    ax[1].set_yscale("log"); ax[1].set_ylabel("Ne (haploid)")
    ax[1].set_xticks(np.arange(len(names)) + .35)
    ax[1].set_xticklabels(names, rotation=20, ha="right", fontsize=9)
    ax[1].legend(P, fontsize=8); ax[1].set_title("fitted Ne per leaf", fontsize=10)
    fig.suptitle(f"Branch drift and Ne  (({P[0]},{P[1]}),{P[2]})", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "drift_and_Ne.png"), dpi=130)
    print(f"[plot] {OUT}/drift_and_Ne.png")


def write_report():
    L = []
    L.append(f"# Three-leaf fit  (({P[0]},{P[1]}),{P[2]})\n")
    L.append(f"engine: {', '.join(ENGINES)}   nodes: {NODES}\n")
    L.append("\n## Parameters\n")
    L.append(f"| model | engine | t_sis | t_root | " + " | ".join(f"Ne_{p}" for p in P)
             + f" | x_{P[0]} | x_{P[1]} | tau | secs |")
    L.append("|" + "---|" * (7 + len(P)))
    for nm, *_ in SPECS:
        for eng in ENGINES:
            f = res.get(f"{nm}|{eng}")
            if not f:
                continue
            L.append(f"| {nm} | {eng} | {f['times'][0]:.0f} | {sum(f['times']):.0f} | "
                     + " | ".join(f"{v:,.0f}" for v in f["Ne"][:3])
                     + f" | {f['x_drift'][0]:.5f} | {f['x_drift'][1]:.5f} "
                     f"| {f['tau']:.2f} | {f['secs']:.0f} |")
    L.append(f"| **OBSERVED** | (w_hat) | | | | | | **{X_OBS[0]:.5f}** | "
             f"**{X_OBS[1]:.5f}** | | |")
    L.append("\n## Likelihood decomposition (per draw)\n")
    L.append("| model | engine | lp_ibd | lp_snp | n_ibd | n_snp | chi2/n ibd | chi2/n snp |")
    L.append("|" + "---|" * 8)
    for nm, *_ in SPECS:
        for eng in ENGINES:
            f = res.get(f"{nm}|{eng}")
            if not f:
                continue
            g = lambda k: (f"{f[k]['mean']:+.1f} ± {f[k]['sd']:.1f}" if k in f else "-")
            n = lambda k: (f"{f[k]['mean']:.0f}" if k in f else "-")
            rr = lambda c, m: (f"{f[c]['mean']/f[m]['mean']:.2f}"
                               if c in f and m in f and f[m]["mean"] else "-")
            L.append(f"| {nm} | {eng} | {g('lp_ibd')} | {g('lp_snp')} | {n('n_ibd_obs')} "
                     f"| {n('n_snp_obs')} | {rr('chi2_ibd','n_ibd_obs')} "
                     f"| {rr('chi2_snp','n_snp_obs')} |")
    L.append("\n## Per-pair IBD fit (Mixed, chi2 per observed bin)\n")
    f = res.get(f"Mixed-Nsm|{ENG0}") or res.get(f"IBD-only-Nsm|{ENG0}")
    if f and "ibd_pred" in f:
        pr = np.asarray(f["ibd_pred"])
        L.append("| pair | bins | chi2/bin | short-bin resid | long-bin resid |")
        L.append("|---|---|---|---|---|")
        for k, (i, j) in enumerate(PAIRS):
            obs, se = IBD_HAT[:, i, j], IBD_SE[:, i, j]
            m = obs > 0
            if m.sum() == 0:
                continue
            r = (obs[m] - pr[m, i, j]) / se[m]
            ns = max(1, m.sum() // 3)
            L.append(f"| {PLBL[k]} | {m.sum()} | {np.mean(r**2):.1f} "
                     f"| {r[:ns].mean():+.1f} | {r[-ns:].mean():+.1f} |")
        L.append("\nA monotone sign change from the short-bin to the long-bin column is the "
                 "signature of Ne changing over time, which constant-per-branch Ne cannot fit.")
    L.append("\n## Figures\n")
    L.append("- `spectrum_fit.png` - observed vs fitted IBD spectrum + residuals, per pair")
    L.append("- `likelihood_split.png` - lp and chi2/n by component, and term counts")
    L.append("- `drift_and_Ne.png` - fitted vs observed SNP branch drift, and Ne per leaf")
    path = os.path.join(OUT, "report.md")
    open(path, "w").write("\n".join(L) + "\n")
    print(f"[report] {path}")


fig_spectrum(); fig_likelihood(); fig_drift(); write_report()
print(f"\n[save] {os.path.join(OUT, 'fits.json')}")
