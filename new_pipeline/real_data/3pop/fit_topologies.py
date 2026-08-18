"""Fit mixed_model_Nsmooth to all 21 three-leaf topologies and rank them by ELBO.

One folder per topology (named by its Newick string) holding spectrum_fit.png,
report.md and fit.json; then comparison/ with the ranking and the selection
probabilities.

------------------------------------------------------------------------------
What is actually being compared
------------------------------------------------------------------------------
Pathfinder is run with psis_resample=False and calculate_lp=True, so every draw
carries both lp__ (log joint) and lp_approx__ (log density of the Gaussian
approximation).  With log w = lp__ - lp_approx__:

    ELBO   = mean(log w)                    lower bound on log Z, gap = KL(q||p)
    logZ   = logsumexp(log w) - log S       importance-sampling estimate of log Z
    ESS    = 1 / sum(w~^2)                  Kish ESS of the self-normalised weights

Report both.  The ELBO is what was asked for, but it is a BOUND, and the bound
gap is a different size for every model -- a model whose posterior the Gaussian
fits badly is penalised for that rather than for fitting the data badly.  logZ
corrects for it to first order.  ESS says whether logZ can be believed at all:
the IS estimate rests on the largest few weights, so a small ESS means the
number is one draw wearing a trenchcoat.  (ESS rather than a Pareto-k tail fit
because it is exact arithmetic rather than an estimator that can itself be
wrong.)

------------------------------------------------------------------------------
The constant correction, and why it is needed here specifically
------------------------------------------------------------------------------
Stan's `~` statements DROP additive constants.  That is harmless when comparing
one model against itself, and it is what lp__ has always been used for here.
It is NOT harmless across models of different dimension, which is exactly what
this script does: the trees have 2 events / 5 nodes, the admixture graphs 4
events / 8 nodes, so they drop DIFFERENT constants.  From the model block:

    times     ~ exponential(0.01)   drops  n_events * log(0.01)
    Ne_raw    ~ std_normal()        drops  n_nodes  * (-0.5*log(2*pi))
    everything else                 identical term count in every model

so  log p_true = lp__ + n_events*log(0.01) - n_nodes*0.5*log(2*pi).

For a tree that is -13.81, for an admixture graph -25.77: the admixture graphs
are over-credited by 11.97 nats before any data is seen.  Applied below by
`lp_const`; the uncorrected numbers are kept as elbo_raw / logz_raw in fit.json
so the correction can be audited.  Within either group it cancels exactly, so
the 18-way admixture ranking is unaffected either way -- it only matters for
tree-vs-admixture, which is the comparison anyone will actually care about.

------------------------------------------------------------------------------
Read the ranking with this in mind
------------------------------------------------------------------------------
An admixture graph does not only add an admixture fraction.  It splits the
admixed leaf's ancestry into two source branches, each with its own Ne, so it
also buys that leaf a PIECEWISE Ne.  We already know these data want exactly
that (the constant-per-branch Ne leaves a monotone tilt in the length-spectrum
residuals).  So a high-scoring admixture graph is evidence for "one leaf's Ne is
not constant" at least as much as for "one leaf is admixed", and the two cannot
be separated on three leaves.  The spectrum residual plots are in each folder
for precisely this reason: look at whether the win came from removing a tilt.
"""
import os, sys, io, json, time, glob, argparse, contextlib, tempfile
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_RD = os.path.dirname(_HERE)
sys.path.insert(0, _RD)
sys.path.insert(0, _HERE)

from cmdstanpy import CmdStanModel                  # noqa: E402
import infer_topology as IT                        # noqa: E402
import fit_dense as FD                             # noqa: E402
import enumerate_3pop as ET                  # noqa: E402

import matplotlib                                  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                    # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True,
                help="Run directory under 3pop/, e.g. gbr_ceu_ibs.  Must match the "
                     "tag given to make_stan_data.sh -- one leaf set per tag, so "
                     "runs never overwrite each other.")
ap.add_argument("--pops", nargs=3, required=True,
                help="MUST be the same order as make_stan_data.sh used: the saved "
                     "matrix rows are indexed by that pop_order.")
ap.add_argument("--prefix", default=None,
                help="Stan input prefix (default 3pop/<tag>/stan_data/stan_<tag>).")
ap.add_argument("--seeds", nargs="+", type=int, default=[1, 7, 13],
                help="Pathfinder is mode-seeking; restarts guard against a bad "
                     "L-BFGS path. The best ELBO over seeds is kept.")
ap.add_argument("--draws", type=int, default=4000)
ap.add_argument("--paths", type=int, default=8)
ap.add_argument("--only", nargs="+", type=int, default=None,
                help="Fit only these topology indices (debugging).")
a = ap.parse_args()

POPS = a.pops
OUT = os.path.join(_HERE, a.tag)
CMP = os.path.join(OUT, "comparison")
os.makedirs(CMP, exist_ok=True)
if a.prefix is None:
    a.prefix = os.path.join(OUT, "stan_data", f"stan_{a.tag}")
# Scratch is per-tag too, so two tags can run concurrently without racing on the
# same CmdStan output directory.
WORK = os.path.join(tempfile.gettempdir(), f"fit_3pop_{a.tag}")

MODEL = CmdStanModel(stan_file=FD._find_model("mixed_model_Nsmooth.stan"))
MODELS = ET.enumerate_all(POPS)
if a.only:
    MODELS = [m for m in MODELS if m["index"] in a.only]

meta = json.load(open(a.prefix + ".json"))
npz = np.load(a.prefix + ".npz")
W_HAT, W_SE = npz["w_hat"], npz["w_se"]
BINS = np.asarray(meta["bins"], dtype=float)
BIN_CTR = BINS.mean(1) if BINS.ndim == 2 else BINS
PAIRS = [(i, j) for i in range(3) for j in range(i, 3)]
PLBL = [f"{POPS[i]}-{POPS[j]}" for i, j in PAIRS]

HALF_LOG_2PI = 0.5 * np.log(2 * np.pi)


def lp_const(n_events, n_nodes):
    """Additive constant Stan's `~` dropped -- see the module docstring."""
    return n_events * np.log(0.01) - n_nodes * HALF_LOG_2PI


def quiet(fn, *args, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kw)


# ---------------------------------------------------------------- fitting
def fit_one(m):
    dem = quiet(ET.make_dem, POPS, m)
    data = quiet(IT.assemble_stan_data, dem, a.prefix + ".npz", a.prefix + ".json",
                 model="mixed", T_max=IT.DEFAULT_T_MAX)
    data.update(FD.smooth_tree_data(dem))
    nodes = list(dem.nodes.keys())
    n_ev, n_ad, n_nd = data["n_events"], data["n_admixture"], data["n_nodes"]
    init = FD.make_init("smooth", False, n_ev, n_ad, n_nd, n_ev + 1)
    C = lp_const(n_ev, n_nd)

    best, runs = None, []
    for seed in a.seeds:
        od = os.path.join(WORK, f"{m['index']:02d}_s{seed}")
        os.makedirs(od, exist_ok=True)
        for f in glob.glob(os.path.join(od, "*")):
            os.remove(f)
        t0 = time.time()
        try:
            fit = MODEL.pathfinder(
                data=data, inits=init, seed=seed, output_dir=od,
                num_paths=a.paths, draws=a.draws, num_single_draws=a.draws // a.paths,
                psis_resample=False, calculate_lp=True, show_console=False)
        except Exception as ex:
            print(f"      seed {seed}: FAILED {str(ex).splitlines()[-1][:70]}", flush=True)
            continue
        cn = list(fit.column_names)
        A = np.asarray(fit.draws()).reshape(-1, len(cn))
        logw = A[:, cn.index("lp__")] - A[:, cn.index("lp_approx__")]
        ok = np.isfinite(logw)
        logw = logw[ok]
        if logw.size < 10:
            print(f"      seed {seed}: only {logw.size} finite draws, skipped", flush=True)
            continue
        S = logw.size
        elbo = float(logw.mean())
        logz = float(np.logaddexp.reduce(logw) - np.log(S))
        wt = np.exp(logw - logw.max()); wt /= wt.sum()
        ess = float(1.0 / np.sum(wt ** 2))
        r = {"seed": seed, "secs": time.time() - t0,
             "elbo_raw": elbo, "logz_raw": logz,
             "elbo": elbo + C, "logz": logz + C,
             "elbo_se": float(logw.std(ddof=1) / np.sqrt(S)),
             "ess": ess, "n_draws": S, "logw": logw}
        sv = {k: np.asarray(v) for k, v in fit.stan_variables().items()}
        r["sv"] = {k: sv[k][ok] if sv[k].shape[0] == ok.size else sv[k] for k in sv}
        print(f"      seed {seed}: ELBO {r['elbo']:+10.1f}  logZ {r['logz']:+10.1f}  "
              f"ESS {ess:7.1f}/{S}  {r['secs']:.0f}s", flush=True)
        runs.append(r)
        # Keep the HIGHEST ELBO across restarts.  Justified, not just optimistic:
        # ELBO <= log Z always, so the largest one found is the tightest bound.
        # But it is a MAXIMUM, so it inherits the seed-to-seed spread -- a model
        # whose L-BFGS paths scatter more gets a higher max for that reason
        # alone.  Every seed's ELBO and log weights are therefore kept, and
        # compare_elbo.py bootstraps over seeds as well as draws so the ranking
        # is judged against that spread rather than in spite of it.
        if best is None or r["elbo"] > best["elbo"]:
            best = r
    if best is None:
        return None

    sv = best.pop("sv")
    for q in runs:
        q.pop("sv", None)
    t = np.asarray(sv["times"]).reshape(len(sv["times"]), -1)
    out = {**{k: v for k, v in best.items() if k != "logw"},
           "index": m["index"], "name": m["name"], "newick": m["newick"],
           "n_admix": m["n_admix"], "admixed": m["admixed"],
           "outside_first_merge_rule": m["outside_first_merge_rule"],
           "nodes": nodes, "n_events": n_ev, "n_nodes": n_nd, "lp_const": float(C),
           "n_seeds": len(runs),
           "elbo_by_seed": {str(q["seed"]): q["elbo"] for q in runs},
           "elbo_seed_min": float(min(q["elbo"] for q in runs)),
           "elbo_seed_max": float(max(q["elbo"] for q in runs)),
           "logw_by_seed": {str(q["seed"]): q["logw"].tolist() for q in runs},
           "events": [dict(e) for e in dem.ordered_events],
           "times": t.mean(0).tolist(),
           "times_sd": t.std(0).tolist(),
           "cum_times": np.cumsum(t, axis=1).mean(0).tolist(),
           "Ne": np.asarray(sv["Ne"]).mean(0).tolist(),
           "Ne_sd_log": np.log(np.asarray(sv["Ne"])).std(0).tolist(),
           "tau": float(np.mean(sv["tau"]))}
    if n_ad:
        f = np.asarray(sv["admixture_fractions"]).reshape(-1, n_ad)
        out["admixture_fractions"] = f.mean(0).tolist()
        out["admixture_fractions_sd"] = f.std(0).tolist()
    for k in ("lp_ibd", "lp_snp", "chi2_ibd", "chi2_snp", "n_ibd_obs", "n_snp_obs"):
        if k in sv:
            v = np.asarray(sv[k]).ravel().astype(float)
            out[k] = float(v.mean())
    if "ibd_fraction" in sv:
        out["ibd_pred"] = np.asarray(sv["ibd_fraction"]).mean(0).tolist()
    if "W_centered" in sv:
        out["W_pred"] = np.asarray(sv["W_centered"]).mean(0).tolist()
    out["logw"] = best["logw"].tolist()
    out["ibd_hat"] = np.asarray(data["ibd_hat"]).tolist()
    out["ibd_se"] = np.asarray(data["ibd_se"]).tolist()
    out["ibd_count"] = np.asarray(data["ibd_count"]).tolist()
    ns = list(data["n_samples"])
    out["n_pairs"] = [[(ns[p] * (ns[p] - 1) / 2.0 if p == q else ns[p] * ns[q])
                       for q in range(len(ns))] for p in range(len(ns))]
    return out


# ---------------------------------------------------------------- per-topology output
def poisson_interval(k, level=0.6827):
    """Exact (Garwood) Poisson interval for a count k, in units of k.

    Returned as multiplicative factors (lo/k, hi/k) so they can scale a rate.
    The plotted quantity is a rate estimated from a count, so its uncertainty is
    the count's.  Two reasons this replaces the +-1 jackknife SE that used to be
    drawn:

      * the jackknife SE collapses to se/mean = sqrt(2/k) in the sparse tail (see
        the model block), so at k = 1 it was 1.414*mean and at k = 2 exactly mean
        -- the lower arm landed at or below ZERO, which has no position on a log
        axis.  That, not anything in the data, is why the lower whiskers ran off
        the bottom of the panel.
      * a symmetric +-se interval is the wrong shape for a count anyway.  The
        Poisson interval is asymmetric the OTHER way (upper arm longer on a log
        axis), which is the honest picture: one observed segment is consistent
        with a rate several times higher far more easily than with one near zero.

    k = 0 gives (0, upper) -- an upper limit, no point estimate.
    """
    from scipy.stats import chi2 as _chi2
    al = 1.0 - level
    k = np.asarray(k, float)
    lo = np.where(k > 0, 0.5 * _chi2.ppf(al / 2, np.maximum(2 * k, 1e-9)), 0.0)
    hi = 0.5 * _chi2.ppf(1 - al / 2, 2 * k + 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(k > 0, lo / k, 0.0), np.where(k > 0, hi / k, np.inf)


def fig_spectrum(r, path):
    """Observed vs fitted IBD length spectrum for all six pairs, plus residuals.

    Top row: rate per pair per cM on a log axis, error bars = exact Poisson.
    The y-range is FIXED and shared across panels -- from 10x the largest value
    anywhere in the figure down to 1e-3 x the smallest -- so the six pairs are
    directly comparable and the model curve stays visible where it falls below
    the last observed segment instead of being clipped at the panel floor.

    Bottom row: Pearson residuals of the Poisson fit, (k - lambda)/sqrt(lambda).
    A MONOTONE tilt across segment length is a leaf whose Ne changed, which an
    admixture graph can absorb through its two source branches without any
    admixture being real -- that is what to look for before believing a win.
    """
    obs = np.asarray(r["ibd_hat"])
    pred = np.asarray(r["ibd_pred"])
    cnt = np.asarray(r["ibd_count"], float)
    npair = np.asarray(r["n_pairs"], float)

    # y-range is set by the DATA, not the model curve.  Top = 10x the largest
    # observed rate.  Bottom = 1e-3 x the smallest rate a data point could take,
    # i.e. the rate contributed by ONE segment in the sparsest pair.  Letting the
    # prediction set the floor would be useless: a tree that cannot make
    # cross-population IBD sends its curve to 1e-36 and squashes everything real
    # into the top 3% of the panel.  The curve running off the bottom is itself
    # the finding, so it is allowed to.
    unit = np.divide(obs, cnt, out=np.full_like(obs, np.inf), where=cnt > 0)
    ylo = unit[np.isfinite(unit)].min() * 1e-3
    yhi = obs.max() * 10.0

    fig, ax = plt.subplots(2, 6, figsize=(23, 7.0), sharex=True)
    for kk, (i, j) in enumerate(PAIRS):
        o, p, c = obs[:, i, j], pred[:, i, j], cnt[:, i, j]
        mk = c > 0
        flo, fhi = poisson_interval(c)
        ax[0, kk].errorbar(BIN_CTR[mk], o[mk],
                           yerr=[o[mk] - o[mk] * flo[mk], o[mk] * fhi[mk] - o[mk]],
                           fmt="o", ms=3, lw=.8, color="k", label="observed", zorder=3)
        ax[0, kk].plot(BIN_CTR, p, lw=1.7, color="#d62728", label="Mixed-Nsm")
        # bins with no segments: draw the Poisson upper limit as a downward arrow
        # bins with no segments carry an upper limit, not a point: Poisson 95%
        # one-sided UL for k = 0 is 3.0 segments, drawn at 3 x the one-segment rate.
        emp = ~mk
        if emp.any() and mk.any():
            u = (o[mk] / c[mk]).mean()
            ax[0, kk].plot(BIN_CTR[emp], np.full(emp.sum(), 2.996 * u),
                           "v", ms=3, color="0.6", label="0 segs (95% UL)")
        ax[0, kk].set_yscale("log")
        ax[0, kk].set_ylim(ylo, yhi)
        ax[0, kk].set_title(f"{PLBL[kk]}   n={c.sum():,.0f} segs", fontsize=10)

        # lambda on the model's own scale: pred is a rate per pair per cM, so
        # lambda = pred * cm * n_pairs.  cm is recovered from the observed side --
        # sum(count) / sum(rate * n_pairs) is exactly cm -- so the residual row
        # needs no extra data plumbed through.
        denom = (o * npair[i, j]).sum()
        cm_hat = (c.sum() / denom) if denom > 0 else np.nan
        lam = p * npair[i, j] * cm_hat
        res = np.where(lam > 0, (c - lam) / np.sqrt(np.maximum(lam, 1e-12)), np.nan)
        ax[1, kk].plot(BIN_CTR, res, ".-", ms=3, lw=1, color="#d62728")
        ax[1, kk].axhline(0, color="k", lw=.8)
        for y in (-2, 2):
            ax[1, kk].axhline(y, color="grey", lw=.6, ls=":")
        ax[1, kk].set_xlabel("segment length (cM)")
        chi = np.nanmean(res ** 2)
        ax[1, kk].text(.97, .04, f"chi2/n = {chi:.1f}", transform=ax[1, kk].transAxes,
                       ha="right", fontsize=9,
                       color=("firebrick" if chi > 8 else "black"))
    ax[0, 0].set_ylabel("IBD rate (per pair per cM)")
    ax[1, 0].set_ylabel("Pearson resid. (k-lam)/sqrt(lam)")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle(f"[{r['index']:02d}] {r['newick']}    ELBO {r['elbo']:+.1f}   "
                 f"logZ {r['logz']:+.1f}   ESS {r['ess']:.0f}", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_report(r, path):
    L = [f"# `{r['newick']}`", ""]
    L.append(f"Model {r['index']:02d} of 21 | "
             f"{'tree (no admixture)' if not r['n_admix'] else 'admixed leaf: **' + r['admixed'] + '**'} | "
             f"{r['n_events']} events, {r['n_nodes']} nodes")
    if r["outside_first_merge_rule"]:
        L.append("")
        L.append("> Graph where the two NON-admixed leaves merge first. Excluded by the "
                 "'first merge must involve an admixture branch' rule; included here "
                 "because it is the standard local-clade + deep-source scenario.")
    L += ["", "## Model score", "",
          "| quantity | value | |", "|---|---|---|",
          f"| ELBO | {r['elbo']:+.2f} | +- {r['elbo_se']:.2f} (MC) |",
          f"| logZ (importance sampling) | {r['logz']:+.2f} | |",
          f"| ESS of the IS weights | {r['ess']:.1f} / {r['n_draws']} | "
          f"{'ok' if r['ess'] >= 50 else 'LOW -- treat logZ with suspicion'} |",
          f"| Stan dropped-constant correction | {r['lp_const']:+.2f} | already applied |",
          f"| seed kept / runtime | {r['seed']} | {r['secs']:.0f} s |", ""]

    L += ["## Events (in temporal order, most recent first)", "",
          "| # | type | detail | time (gen) | cumulative |", "|---|---|---|---|---|"]
    for k, ev in enumerate(r["events"]):
        if ev["type"] == "MERGE":
            det = f"{ev['children'][0]} + {ev['children'][1]} -> {ev['parent']}"
        else:
            det = f"{ev['child']} -> {ev['parents'][0]} + {ev['parents'][1]}"
        L.append(f"| {k+1} | {ev['type']} | {det} | "
                 f"{r['times'][k]:,.1f} +- {r['times_sd'][k]:,.1f} | "
                 f"{r['cum_times'][k]:,.1f} |")
    L.append("")

    if r["n_admix"]:
        f = r["admixture_fractions"][0]
        L += ["## Admixture fraction", "",
              f"**f = {f:.3f} +- {r['admixture_fractions_sd'][0]:.3f}** "
              f"(fraction from `{r['events'][0]['parents'][0]}`; "
              f"{1-f:.3f} from `{r['events'][0]['parents'][1]}`)", ""]
        if min(f, 1 - f) < 0.05:
            L.append("Collapsed to a tree: one source carries <5% of the ancestry, so "
                     "this graph is behaving as its no-admixture special case.")
            L.append("")

    L += ["## Effective sizes (haploid)", "",
          "| node | Ne | sd of log Ne |", "|---|---|---|"]
    for nm, ne, sd in zip(r["nodes"], r["Ne"], r["Ne_sd_log"]):
        L.append(f"| `{nm}` | {ne:,.0f} | {sd:.2f} |")
    L += ["", f"log-Ne random-walk step scale tau = {r['tau']:.3f}", ""]

    L += ["## Fit quality by component", "",
          "| component | log-likelihood | n terms | chi2/n |", "|---|---|---|---|"]
    for c, lp, ch, n in (("IBD", "lp_ibd", "chi2_ibd", "n_ibd_obs"),
                         ("SNP", "lp_snp", "chi2_snp", "n_snp_obs")):
        if lp in r:
            L.append(f"| {c} | {r[lp]:+,.1f} | {r[n]:,.0f} | {r[ch]/max(r[n],1):.2f} |")
    L += ["", "chi2/n is the mean squared standardized residual: 1.0 = fits inside its "
          "own SEs. Do NOT compare the two log-likelihood LEVELS against each other -- "
          "different term counts and different units.", ""]

    if "W_pred" in r:
        Wp = np.asarray(r["W_pred"])
        L += ["## SNP covariance residuals `(w_hat - W_pred)/w_se`", "",
              "| | " + " | ".join(POPS) + " |", "|---|" + "---|" * 3]
        for i in range(3):
            row = [f"{(W_HAT[i,j]-Wp[i,j])/W_SE[i,j]:+.2f}" for j in range(3)]
            L.append(f"| **{POPS[i]}** | " + " | ".join(row) + " |")
        L.append("")
    L.append("![spectrum](spectrum_fit.png)")
    open(path, "w").write("\n".join(L) + "\n")


# ---------------------------------------------------------------- run
print(f"[3pop] {len(MODELS)} topologies on {POPS}\n")
results = []
for m in MODELS:
    print(f"[{m['index']:02d}/21] {m['newick']}", flush=True)
    r = fit_one(m)
    if r is None:
        print("      ALL SEEDS FAILED", flush=True)
        continue
    d = os.path.join(OUT, m["name"])
    os.makedirs(d, exist_ok=True)
    if "ibd_pred" in r:
        fig_spectrum(r, os.path.join(d, "spectrum_fit.png"))
    write_report(r, os.path.join(d, "report.md"))
    BIG = ("logw", "logw_by_seed")
    json.dump({k: v for k, v in r.items() if k not in BIG},
              open(os.path.join(d, "fit.json"), "w"), indent=2)
    results.append(r)

json.dump({"pops": POPS, "results": [{k: v for k, v in r.items()
                                      if k not in ("logw", "logw_by_seed",
                                                   "ibd_hat", "ibd_se")}
                                     for r in results]},
          open(os.path.join(CMP, "all_fits.json"), "w"), indent=2)
np.savez(os.path.join(CMP, "logw.npz"),
         **{f"{r['index']}_s{s}": np.asarray(v)
            for r in results for s, v in r["logw_by_seed"].items()})
print(f"\n[done] {len(results)}/{len(MODELS)} fitted -> now run compare_elbo.py")
