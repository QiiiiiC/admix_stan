"""Fit mixed Nsmooth variants to all 21 three-leaf topologies.

Each topology contains one subdirectory per variant, holding spectrum_fit.png,
report.md and fit.json. Variant-wide fit summaries are written below comparison/.
The original four variants are the default; four optional recent-grid variants
can be selected with --variants.

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
ap.add_argument("--variants", nargs="+", default=None,
                help="Subset of model keys to fit (default: all four).")
ap.add_argument("--only", nargs="+", type=int, default=None,
                help="Fit only these topology indices (debugging).")
a = ap.parse_args()

POPS = a.pops
OUT = os.path.join(_HERE, a.tag)
os.makedirs(os.path.join(OUT, "comparison"), exist_ok=True)
if a.prefix is None:
    a.prefix = os.path.join(OUT, "stan_data", f"stan_{a.tag}")
# Scratch is per-tag too, so two tags can run concurrently without racing on the
# same CmdStan output directory.
WORK = os.path.join(tempfile.gettempdir(), f"fit_3pop_{a.tag}")

SPECS = {
    "poisson_shared_ne": {
        "stan": "mixed_model_Nsmooth_poisson.stan",
        "likelihood": "poisson", "separate_ne": False,
        "label": "Poisson, shared Ne",
    },
    "normal_shared_ne": {
        "stan": "mixed_model_Nsmooth_normal.stan",
        "likelihood": "normal", "separate_ne": False,
        "label": "Normal, shared Ne",
    },
    "poisson_separate_ne": {
        "stan": "mixed_model_Nsmooth_poisson_separate_ne.stan",
        "likelihood": "poisson", "separate_ne": True,
        "label": "Poisson, separate IBD/SNP Ne",
    },
    "normal_separate_ne": {
        "stan": "mixed_model_Nsmooth_normal_separate_ne.stan",
        "likelihood": "normal", "separate_ne": True, "grid": False,
        "label": "Normal, separate IBD/SNP Ne",
    },
    "poisson_grid_shared_ne": {
        "stan": "mixed_model_Ngrid_poisson.stan",
        "likelihood": "poisson", "separate_ne": False, "grid": True,
        "label": "Poisson, shared Ne, recent grid",
    },
    "normal_grid_shared_ne": {
        "stan": "mixed_model_Ngrid_normal.stan",
        "likelihood": "normal", "separate_ne": False, "grid": True,
        "label": "Normal, shared Ne, recent grid",
    },
    "poisson_grid_separate_ne": {
        "stan": "mixed_model_Ngrid_poisson_separate_ne.stan",
        "likelihood": "poisson", "separate_ne": True, "grid": True,
        "label": "Poisson, separate IBD/SNP Ne, recent grid",
    },
    "normal_grid_separate_ne": {
        "stan": "mixed_model_Ngrid_normal_separate_ne.stan",
        "likelihood": "normal", "separate_ne": True, "grid": True,
        "label": "Normal, separate IBD/SNP Ne, recent grid",
    },
}

DEFAULT_VARIANTS = [
    "poisson_shared_ne", "normal_shared_ne",
    "poisson_separate_ne", "normal_separate_ne",
]
for spec in SPECS.values():
    spec.setdefault("grid", False)


def find_model(filename):
    path = os.path.join(_RD, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def smooth_tree_data(dem):
    """Encode the Nsmooth parent/random-walk structure in node order."""
    names = list(dem.nodes.keys())
    idx = {name: i + 1 for i, name in enumerate(names)}
    ne_parent = [0] * len(names)
    ne_admix_idx = [0] * len(names)
    ne_start_event = [0] * len(names)
    admixture_index = 0
    for event_index, event in enumerate(dem.ordered_events, start=1):
        if event["type"] == "MERGE":
            parent = event["parent"]
            child1, child2 = event["children"]
            ne_start_event[idx[parent] - 1] = event_index
            ne_parent[idx[child1] - 1] = idx[parent]
            ne_parent[idx[child2] - 1] = idx[parent]
        else:
            admixture_index += 1
            child = event["child"]
            source1, source2 = event["parents"]
            ne_start_event[idx[source1] - 1] = event_index
            ne_start_event[idx[source2] - 1] = event_index
            ne_admix_idx[idx[child] - 1] = admixture_index
    return {
        "ne_parent": ne_parent,
        "ne_admix_idx": ne_admix_idx,
        "ne_start_event": ne_start_event,
    }


def shared_ne_init(n_events, n_admixture, n_nodes, n_leaves, grid=False):
    init = {
        "times": [100.0] * n_events,
        "admixture_fractions": [0.5] * n_admixture,
        "mu_log": float(np.log(15000)),
        "sigma_log": 0.3,
        "tau": 0.3,
        "Ne_raw": [0.0] * n_nodes,
    }
    if grid:
        init["Ne_recent_raw"] = [[0.0, 0.0] for _ in range(n_leaves)]
    return init


def shared_grid_init_from_base(m, data, n_events, n_admixture, n_nodes,
                               n_leaves, likelihood):
    """Place a shared-grid model on its fitted non-grid submodel."""
    base_variant = f"{likelihood}_shared_ne"
    path = os.path.join(OUT, m["name"], base_variant, "fit.json")
    if not os.path.exists(path):
        print(f"      [init] {base_variant} fit unavailable; using prior center", flush=True)
        return shared_ne_init(n_events, n_admixture, n_nodes, n_leaves, grid=True)

    base = json.load(open(path))
    event_durations = np.maximum(np.asarray(base["times"], float), 1.001)
    event_durations[0] = max(event_durations[0], 11.001)
    sampled_times = event_durations.copy()
    sampled_times[0] -= 10.0
    fractions = np.clip(
        np.asarray(base.get("admixture_fractions", [0.5] * n_admixture), float),
        1e-3, 1.0 - 1e-3,
    )
    log_ne = np.log(np.asarray(base["Ne"], float))
    tau = max(float(base["tau"]), 0.05)
    sigma = 0.3

    ne_parent = np.asarray(data["ne_parent"], int)
    ne_admix_idx = np.asarray(data["ne_admix_idx"], int)
    start_event = np.asarray(data["ne_start_event"], int)
    cumulative_times = np.cumsum(event_durations)
    node_t = np.where(start_event > 0,
                      cumulative_times[np.maximum(start_event - 1, 0)], 10.0)
    roots = (ne_parent == 0) & (ne_admix_idx == 0)
    mu = float(log_ne[roots].mean())
    raw = np.zeros(n_nodes)
    admixture_map = np.asarray(data["admixture_map"], int)
    for node in range(n_nodes - 1, -1, -1):
        if roots[node]:
            raw[node] = (log_ne[node] - mu) / sigma
        elif ne_admix_idx[node] > 0:
            admixture = ne_admix_idx[node] - 1
            source1 = admixture_map[admixture, 2]
            source2 = admixture_map[admixture, 3]
            center = (fractions[admixture] * log_ne[source1]
                      + (1.0 - fractions[admixture]) * log_ne[source2])
            dt = max(node_t[source1] - node_t[node], 1e-7)
            raw[node] = ((log_ne[node] - center)
                         / (tau * np.sqrt(dt / 100.0 + 1e-9)))
        else:
            parent = ne_parent[node] - 1
            dt = max(node_t[parent] - node_t[node], 1e-7)
            raw[node] = ((log_ne[node] - log_ne[parent])
                         / (tau * np.sqrt(dt / 100.0 + 1e-9)))
    if np.max(np.abs(raw)) > 6.0:
        print(f"      [init] clipping extreme Ne_raw coordinate "
              f"({np.max(np.abs(raw)):.1f} -> 6)", flush=True)
        raw = np.clip(raw, -6.0, 6.0)

    print(f"      [init] warm start from {base_variant} (ELBO {base['elbo']:+.1f})",
          flush=True)
    return {
        "times": sampled_times.tolist(),
        "admixture_fractions": fractions.tolist(),
        "mu_log": mu,
        "sigma_log": sigma,
        "tau": tau,
        "Ne_raw": raw.tolist(),
        "Ne_recent_raw": [[0.0, 0.0] for _ in range(n_leaves)],
    }


def separate_ne_init_from_shared(m, data, n_events, n_admixture, n_nodes,
                                 n_leaves, likelihood, grid=False):
    """Start the larger model on the fitted shared-Ne manifold."""
    shared_variant = f"{likelihood}_{'grid_' if grid else ''}shared_ne"
    path = os.path.join(OUT, m["name"], shared_variant, "fit.json")
    if not os.path.exists(path):
        print(f"      [init] {shared_variant} fit unavailable; using prior center", flush=True)
        init = {
            "times": [100.0] * n_events,
            "admixture_fractions": [0.5] * n_admixture,
            "mu_log_ibd": float(np.log(15000)), "sigma_log_ibd": 0.3,
            "tau_ibd": 0.3, "Ne_raw_ibd": [0.0] * n_nodes,
            "mu_log_snp": float(np.log(15000)), "sigma_log_snp": 0.3,
            "tau_snp": 0.3, "Ne_raw_snp": [0.0] * n_nodes,
        }
        if grid:
            for component in ("ibd", "snp"):
                init[f"Ne_recent_raw_{component}"] = [
                    [0.0, 0.0] for _ in range(n_leaves)]
        return init

    shared = json.load(open(path))
    event_durations = np.maximum(np.asarray(shared["times"], float), 1.001)
    times = event_durations.copy()
    if grid:
        times[0] = max(event_durations[0] - 10.0, 1.001)
    fractions = np.clip(
        np.asarray(shared.get("admixture_fractions", [0.5] * n_admixture), float),
        1e-3, 1.0 - 1e-3,
    )
    log_ne = np.log(np.asarray(shared["Ne"], float))
    tau = max(float(shared["tau"]), 0.05)
    sigma = 0.3

    ne_parent = np.asarray(data["ne_parent"], int)
    ne_admix_idx = np.asarray(data["ne_admix_idx"], int)
    start_event = np.asarray(data["ne_start_event"], int)
    cumulative_times = np.cumsum(event_durations)
    leaf_start = 10.0 if grid else 0.0
    node_t = np.where(start_event > 0,
                      cumulative_times[np.maximum(start_event - 1, 0)], leaf_start)
    roots = (ne_parent == 0) & (ne_admix_idx == 0)
    mu = float(log_ne[roots].mean())
    raw = np.zeros(n_nodes)
    admixture_map = np.asarray(data["admixture_map"], int)
    for node in range(n_nodes - 1, -1, -1):
        if roots[node]:
            raw[node] = (log_ne[node] - mu) / sigma
        elif ne_admix_idx[node] > 0:
            admixture = ne_admix_idx[node] - 1
            source1 = admixture_map[admixture, 2]
            source2 = admixture_map[admixture, 3]
            center = (fractions[admixture] * log_ne[source1]
                      + (1.0 - fractions[admixture]) * log_ne[source2])
            dt = max(node_t[source1] - node_t[node], 1e-7)
            raw[node] = (log_ne[node] - center) / (tau * np.sqrt(dt / 100.0 + 1e-9))
        else:
            parent = ne_parent[node] - 1
            dt = max(node_t[parent] - node_t[node], 1e-7)
            raw[node] = ((log_ne[node] - log_ne[parent])
                         / (tau * np.sqrt(dt / 100.0 + 1e-9)))

    # Posterior means of transformed Ne values need not map back to a plausible
    # mean noncentered coordinate.  This especially matters for ancestry paths
    # whose fitted admixture weight is essentially zero: their Ne can wander to
    # extreme values without changing the likelihood.  Keep the warm start in a
    # finite, prior-plausible region while retaining the shared fit elsewhere.
    if np.max(np.abs(raw)) > 6.0:
        print(f"      [init] clipping extreme Ne_raw coordinate "
              f"({np.max(np.abs(raw)):.1f} -> 6)", flush=True)
        raw = np.clip(raw, -6.0, 6.0)

    init = {"times": times.tolist(), "admixture_fractions": fractions.tolist()}
    for component in ("ibd", "snp"):
        init[f"mu_log_{component}"] = mu
        init[f"sigma_log_{component}"] = sigma
        init[f"tau_{component}"] = tau
        init[f"Ne_raw_{component}"] = raw.tolist()
        if grid:
            recent = np.asarray(shared["Ne_recent"], float)
            recent_raw = np.zeros((n_leaves, 2))
            scale = tau * np.sqrt(5.0 / 100.0)
            recent_raw[:, 1] = (np.log(recent[:, 1]) - log_ne[:n_leaves]) / scale
            recent_raw[:, 0] = (np.log(recent[:, 0]) - np.log(recent[:, 1])) / scale
            init[f"Ne_recent_raw_{component}"] = np.clip(
                recent_raw, -6.0, 6.0).tolist()
    print(f"      [init] warm start from {shared_variant} (ELBO {shared['elbo']:+.1f})", flush=True)
    return init


variant_keys = DEFAULT_VARIANTS if a.variants is None else a.variants
unknown = sorted(set(variant_keys) - set(SPECS))
if unknown:
    ap.error(f"unknown variants: {unknown}; choose from {list(SPECS)}")
for key in variant_keys:
    SPECS[key]["key"] = key
    SPECS[key]["model"] = CmdStanModel(stan_file=find_model(SPECS[key]["stan"]))

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


def lp_const(n_events, n_nodes, n_leaves, separate_ne, grid=False):
    """Restore dropped prior constants for topology and Ne-model comparison."""
    n_traj = 2 if separate_ne else 1
    c_mu = -np.log(0.53) - HALF_LOG_2PI
    # sigma_log and tau are constrained positive and intended as half-Normals.
    c_half = np.log(2.0) - np.log(0.3) - HALF_LOG_2PI
    n_raw = n_nodes + (2 * n_leaves if grid else 0)
    c_traj = c_mu + 2 * c_half - n_raw * HALF_LOG_2PI
    return n_events * np.log(0.01) + n_traj * c_traj


def quiet(fn, *args, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kw)


# ---------------------------------------------------------------- fitting
def fit_one(m, spec):
    dem = quiet(ET.make_dem, POPS, m)
    data = quiet(IT.assemble_stan_data, dem, a.prefix + ".npz", a.prefix + ".json",
                 model="mixed", T_max=IT.DEFAULT_T_MAX)
    data.update(smooth_tree_data(dem))
    nodes = list(dem.nodes.keys())
    n_ev, n_ad, n_nd = data["n_events"], data["n_admixture"], data["n_nodes"]
    n_leaves = data["n_leaves"]
    if spec["separate_ne"]:
        init = separate_ne_init_from_shared(
            m, data, n_ev, n_ad, n_nd, n_leaves, spec["likelihood"], spec["grid"])
        cold_init = None
    else:
        if spec["grid"]:
            init = shared_grid_init_from_base(
                m, data, n_ev, n_ad, n_nd, n_leaves, spec["likelihood"])
            cold_init = shared_ne_init(n_ev, n_ad, n_nd, n_leaves, grid=True)
        else:
            init = shared_ne_init(n_ev, n_ad, n_nd, n_leaves)
            cold_init = None
    C = lp_const(n_ev, n_nd, n_leaves, spec["separate_ne"], spec["grid"])

    best, runs = None, []
    for seed_index, seed in enumerate(a.seeds):
        od = os.path.join(WORK, spec["key"], f"{m['index']:02d}_s{seed}")
        os.makedirs(od, exist_ok=True)
        for f in glob.glob(os.path.join(od, "*")):
            os.remove(f)
        t0 = time.time()
        run_init = (cold_init if cold_init is not None and
                    seed_index == len(a.seeds) - 1 else init)
        if run_init is cold_init and cold_init is not None:
            print(f"      seed {seed}: prior-center grid start", flush=True)
        try:
            fit = spec["model"].pathfinder(
                data=data, inits=run_init, seed=seed, output_dir=od,
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
    time_variable = "event_durations" if spec["grid"] else "times"
    t = np.asarray(sv[time_variable]).reshape(len(sv[time_variable]), -1)
    out = {**{k: v for k, v in best.items() if k != "logw"},
           "index": m["index"], "name": m["name"], "newick": m["newick"],
           "variant": spec["key"], "variant_label": spec["label"],
           "likelihood": spec["likelihood"], "separate_ne": spec["separate_ne"],
           "recent_grid": spec["grid"],
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
           "cum_times": np.cumsum(t, axis=1).mean(0).tolist()}
    if spec["separate_ne"]:
        for component in ("ibd", "snp"):
            ne = np.asarray(sv[f"Ne_{component}"])
            out[f"Ne_{component}"] = ne.mean(0).tolist()
            out[f"Ne_{component}_sd_log"] = np.log(ne).std(0).tolist()
            out[f"mu_log_{component}"] = float(np.mean(sv[f"mu_log_{component}"]))
            out[f"sigma_log_{component}"] = float(np.mean(sv[f"sigma_log_{component}"]))
            out[f"tau_{component}"] = float(np.mean(sv[f"tau_{component}"]))
            out[f"Ne_raw_{component}"] = np.asarray(
                sv[f"Ne_raw_{component}"]).mean(0).tolist()
            if spec["grid"]:
                recent = np.asarray(sv[f"Ne_recent_{component}"])
                out[f"Ne_recent_{component}"] = recent.mean(0).tolist()
                out[f"Ne_recent_{component}_sd_log"] = np.log(recent).std(0).tolist()
                out[f"Ne_recent_raw_{component}"] = np.asarray(
                    sv[f"Ne_recent_raw_{component}"]).mean(0).tolist()
    else:
        ne = np.asarray(sv["Ne"])
        out["Ne"] = ne.mean(0).tolist()
        out["Ne_sd_log"] = np.log(ne).std(0).tolist()
        out["mu_log"] = float(np.mean(sv["mu_log"]))
        out["sigma_log"] = float(np.mean(sv["sigma_log"]))
        out["tau"] = float(np.mean(sv["tau"]))
        out["Ne_raw"] = np.asarray(sv["Ne_raw"]).mean(0).tolist()
        if spec["grid"]:
            recent = np.asarray(sv["Ne_recent"])
            out["Ne_recent"] = recent.mean(0).tolist()
            out["Ne_recent_sd_log"] = np.log(recent).std(0).tolist()
            out["Ne_recent_raw"] = np.asarray(sv["Ne_recent_raw"]).mean(0).tolist()
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
    if "ibd_number" in sv:
        out["ibd_number_pred"] = np.asarray(sv["ibd_number"]).mean(0).tolist()
    if "ibd_theory_se" in sv:
        out["ibd_theory_se_pred"] = np.asarray(sv["ibd_theory_se"]).mean(0).tolist()
    if "W_centered" in sv:
        out["W_pred"] = np.asarray(sv["W_centered"]).mean(0).tolist()
    out["logw"] = best["logw"].tolist()
    out["ibd_hat"] = np.asarray(data["ibd_hat"]).tolist()
    out["ibd_se"] = np.asarray(data["ibd_se"]).tolist()
    out["ibd_count"] = np.asarray(data["ibd_count"]).tolist()
    out["cm"] = float(data["cm"])
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
    """Likelihood-specific IBD spectrum and standardized residuals."""
    obs = np.asarray(r["ibd_hat"])
    pred = np.asarray(r["ibd_pred"])
    cnt = np.asarray(r["ibd_count"], float)
    npair = np.asarray(r["n_pairs"], float)
    cm = float(r["cm"])

    fig, ax = plt.subplots(2, 6, figsize=(23, 7.0), sharex=True)
    for kk, (i, j) in enumerate(PAIRS):
        o, p, c = obs[:, i, j], pred[:, i, j], cnt[:, i, j]
        if r["likelihood"] == "poisson":
            q_obs = c / (cm * npair[i, j])
            q_pred = np.asarray(r["ibd_number_pred"])[:, i, j]
            flo, fhi = poisson_interval(c)
            mk = c > 0
            ax[0, kk].errorbar(
                BIN_CTR[mk], q_obs[mk],
                yerr=[q_obs[mk] - q_obs[mk] * flo[mk],
                      q_obs[mk] * fhi[mk] - q_obs[mk]],
                fmt="o", ms=3, lw=.8, color="k", label="observed count rate",
                zorder=3,
            )
            emp = ~mk
            if emp.any():
                upper = 2.996 / (cm * npair[i, j])
                ax[0, kk].plot(BIN_CTR[emp], np.full(emp.sum(), upper), "v",
                               ms=3, color="0.6", label="0 count (95% UL)")
            ax[0, kk].plot(BIN_CTR, q_pred, lw=1.7, color="#d62728",
                           label="expected count rate")
            positives = q_obs[mk]
            if positives.size:
                ax[0, kk].set_yscale("log")
                ax[0, kk].set_ylim(positives.min() * 0.1, positives.max() * 10)
            lam = np.maximum(q_pred * cm * npair[i, j], 1e-12)
            res = (c - lam) / np.sqrt(lam)
        else:
            se = np.asarray(r["ibd_theory_se_pred"])[:, i, j]
            ax[0, kk].errorbar(BIN_CTR, o, yerr=se, fmt="o", ms=3, lw=.8,
                               color="k", label="observed fraction")
            ax[0, kk].plot(BIN_CTR, p, lw=1.7, color="#d62728",
                           label="expected fraction")
            ax[0, kk].fill_between(BIN_CTR, np.maximum(p - se, 0), p + se,
                                   color="#d62728", alpha=.12, linewidth=0)
            res = (o - p) / np.maximum(se, 1e-12)

        ax[0, kk].set_title(f"{PLBL[kk]}   n={c.sum():,.0f} segs", fontsize=10)
        ax[1, kk].plot(BIN_CTR, res, ".-", ms=3, lw=1, color="#d62728")
        ax[1, kk].axhline(0, color="k", lw=.8)
        for y in (-2, 2):
            ax[1, kk].axhline(y, color="grey", lw=.6, ls=":")
        ax[1, kk].set_xlabel("segment length (cM)")
        chi = np.nanmean(res ** 2)
        ax[1, kk].text(.97, .04, f"chi2/n = {chi:.1f}", transform=ax[1, kk].transAxes,
                       ha="right", fontsize=9,
                       color=("firebrick" if chi > 8 else "black"))
    ax[0, 0].set_ylabel("segment rate per pair per cM" if r["likelihood"] == "poisson"
                        else "mean IBD fraction")
    ax[1, 0].set_ylabel("Pearson residual" if r["likelihood"] == "poisson"
                        else "Normal residual")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle(f"{r['variant_label']} | [{r['index']:02d}] {r['newick']}    "
                 f"ELBO {r['elbo']:+.1f}   "
                 f"logZ {r['logz']:+.1f}   ESS {r['ess']:.0f}", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_report(r, path):
    L = [f"# `{r['newick']}`", ""]
    L.append(f"**{r['variant_label']}** | topology {r['index']:02d} of 21 | "
             f"{'tree (no admixture)' if not r['n_admix'] else 'admixed leaf: **' + r['admixed'] + '**'} | "
             f"{r['n_events']} events, {r['n_nodes']} nodes")
    if r.get("recent_grid"):
        L += ["", "> Recent Ne grid: 0-5 and 5-10 generations. The first "
              "demographic event is explicitly constrained to occur after "
              "generation 10 (minimum 11 with the current one-generation gap)."]
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

    if r["separate_ne"]:
        for component in ("ibd", "snp"):
            if r.get("recent_grid"):
                L += [f"## Recent effective sizes for {component.upper()} (haploid)", "",
                      "| population | 0-5 gen | 5-10 gen |", "|---|---:|---:|"]
                for pop, values in zip(POPS, r[f"Ne_recent_{component}"]):
                    L.append(f"| `{pop}` | {values[0]:,.0f} | {values[1]:,.0f} |")
                L.append("")
            L += [f"## Effective sizes for {component.upper()} (haploid)", "",
                  "| node | Ne | sd of log Ne |", "|---|---|---|"]
            for nm, ne, sd in zip(r["nodes"], r[f"Ne_{component}"],
                                  r[f"Ne_{component}_sd_log"]):
                L.append(f"| `{nm}` | {ne:,.0f} | {sd:.2f} |")
            L += ["", f"log-Ne random-walk step scale tau_{component} = "
                        f"{r[f'tau_{component}']:.3f}", ""]
    else:
        if r.get("recent_grid"):
            L += ["## Recent effective sizes shared by IBD and SNP (haploid)", "",
                  "| population | 0-5 gen | 5-10 gen |", "|---|---:|---:|"]
            for pop, values in zip(POPS, r["Ne_recent"]):
                L.append(f"| `{pop}` | {values[0]:,.0f} | {values[1]:,.0f} |")
            L.append("")
        L += ["## Effective sizes shared by IBD and SNP (haploid)", "",
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
    ibd_noise = ("Poisson counting variance" if r["likelihood"] == "poisson"
                 else "Palamara model-derived Normal variance")
    L += ["", f"IBD chi2/n uses {ibd_noise}; SNP chi2/n uses its block SEs. "
          "Values near 1 indicate residuals on the modeled noise scale.", ""]

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
print(f"[3pop] {len(MODELS)} topologies x {len(variant_keys)} variants on {POPS}\n")
for variant in variant_keys:
    spec = SPECS[variant]
    print(f"\n=== {spec['label']} ===", flush=True)
    results = []
    for m in MODELS:
        print(f"[{m['index']:02d}/21] {m['newick']}", flush=True)
        r = fit_one(m, spec)
        if r is None:
            print("      ALL SEEDS FAILED", flush=True)
            continue
        d = os.path.join(OUT, m["name"], variant)
        os.makedirs(d, exist_ok=True)
        if "ibd_pred" in r:
            fig_spectrum(r, os.path.join(d, "spectrum_fit.png"))
        write_report(r, os.path.join(d, "report.md"))
        BIG = ("logw", "logw_by_seed")
        with open(os.path.join(d, "fit.json"), "w") as fh:
            json.dump({k: v for k, v in r.items() if k not in BIG}, fh, indent=2)
        results.append(r)

    cmp = os.path.join(OUT, "comparison", variant)
    os.makedirs(cmp, exist_ok=True)
    compact = [{k: v for k, v in r.items()
                if k not in ("logw", "logw_by_seed", "ibd_hat", "ibd_se")}
               for r in results]
    all_fits_path = os.path.join(cmp, "all_fits.json")
    if a.only and os.path.exists(all_fits_path):
        previous = json.load(open(all_fits_path)).get("results", [])
        compact = list({r["index"]: r for r in previous + compact}.values())
        compact.sort(key=lambda r: r["index"])
    with open(all_fits_path, "w") as fh:
        json.dump({"pops": POPS, "variant": variant, "spec": {
            k: v for k, v in spec.items() if k != "model"
        }, "results": compact}, fh, indent=2)
    logw_path = os.path.join(cmp, "logw.npz")
    weights = {}
    if a.only and os.path.exists(logw_path):
        old_weights = np.load(logw_path)
        weights.update({k: old_weights[k] for k in old_weights.files})
    weights.update({f"{r['index']}_s{s}": np.asarray(v)
                    for r in results for s, v in r["logw_by_seed"].items()})
    np.savez(logw_path, **weights)
    print(f"[done] {len(results)}/{len(MODELS)} fitted for {variant}")
