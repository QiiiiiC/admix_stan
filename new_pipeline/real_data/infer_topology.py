"""
Minimal real-data pipeline:  topology + SNP/IBD data  ->  parameters along the
topology.

Given
  * a candidate ``DemographicTopology`` over the 1000G superpopulations
    (see ``topologies.py``), and
  * the Stan data product written by ``generate_stan_data.py``
    (``<prefix>.npz`` + ``<prefix>.json`` = w_hat/w_se and ibd_mean/ibd_var),

this assembles a cmdstanpy data dict (reusing ``assemble_stan_data``), fits the
chosen Stan model with Pathfinder, and reports the posterior-mean estimate of
every parameter mapped back onto the topology:

  * one event time per MERGE / ADMIXTURE event,
  * one admixture fraction per ADMIXTURE event,
  * effective sizes (per-node for the Nvarying models, a single shared
    ``effective_N`` for the Nfixed models).

Outputs (``--out`` prefix, default next to the data):
  <prefix>_estimates.json  : labelled parameter table (means + SDs)
  <prefix>_topology.png    : the topology drawn with the fitted times/fractions

Defaults match the choices for this project: model=mixed, N=Nvarying, Pathfinder.

Examples
--------
PY=/opt/miniconda3/envs/genetics_env/bin/python

  # 1) consume an existing data product
  $PY infer_topology.py --data real_stan_data --topology amr_admixture

  # 2) end-to-end: build the SNP+IBD product first (needs an IBD caller), then fit
  $PY infer_topology.py --generate --ibd-method hapibd --chrom 22 \
       --data real_chr22 --topology amr_admixture

The IBD caller for ``--generate`` is whatever ``generate_stan_data.py`` supports
(hap-IBD needs the JAR on PATH or HAPIBD_CMD set; FastSMC needs the asmc package).
"""

import os
import sys
import json
import argparse
import subprocess

import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_THIS_DIR)
for _p in (_THIS_DIR, _PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from cmdstanpy import CmdStanModel  # noqa: E402

from generate_stan_data import assemble_stan_data  # noqa: E402
from topologies import load_topology, DEFAULT_POP_ORDER  # noqa: E402

_MODELS = os.path.join(_PARENT, "models")

# Default upper bound for the final open-ended IBD epoch (generations), matching
# the simulation run files.
DEFAULT_T_MAX = 100000.0


# ----------------------------------------------------------------------
# Pathfinder init (mirrors the run files so Pathfinder starts near sane values)
# ----------------------------------------------------------------------
def build_init(dem, model, varying):
    """Init dict for Pathfinder, matching the simulation run files."""
    n_events = len(dem.ordered_events)
    n_nodes = len(dem.nodes)
    n_admix = dem.n_admix

    init = {"times": [100.0] * n_events}
    if n_admix > 0:
        init["admixture_fractions"] = [0.5] * n_admix
    if varying:
        # Nvarying: non-centered hierarchical log-normal Ne.
        #   Ne[a] = exp(mu_log + sigma_log * Ne_raw[a])
        # Ne_raw is a STANDARD-NORMAL deviate (NOT Ne itself!), so it must be
        # O(1); the Ne level is set by mu_log.  Ne_raw=0 => Ne=exp(mu_log)=15000.
        # (Initialising Ne_raw to e.g. 10000 makes Ne=exp(3000)=inf and every
        # Pathfinder path dies with a non-finite gradient.)
        init["mu_log"] = float(np.log(15000.0))
        init["sigma_log"] = 0.3
        init["Ne_raw"] = [0.0] * n_nodes
    else:
        # Nfixed: single shared effective size.
        init["effective_N"] = 15000.0
    return init


def stan_file_for(model, varying):
    suffix = "Nvarying" if varying else "Nfixed"
    return os.path.join(_MODELS, f"{model}_model_{suffix}.stan")


# ----------------------------------------------------------------------
# Map fitted posterior means onto the topology
# ----------------------------------------------------------------------
def extract_estimates(dem, sv, varying):
    """Turn fit.stan_variables() into a labelled estimate dict.

    Conventions (identical to relative_error.build_param_spec):
      * cumulative_times[e]      -> time of ordered event e (MERGE: parent
                                    start; ADMIXTURE: child end),
      * admixture_fractions[i]   -> fraction to the FIRST parent of the i-th
                                    admixture (chronological),
      * Ne[k] (haploid)          -> node list(dem.nodes)[k]  (Nvarying), or the
                                    single effective_N for every node (Nfixed).
    """
    def col(name, idx):
        arr = np.asarray(sv[name], dtype=float)
        v = arr.mean() if arr.ndim == 1 else arr[:, idx].mean()
        s = arr.std() if arr.ndim == 1 else arr[:, idx].std()
        return float(v), float(s)

    node_names = list(dem.nodes.keys())
    cum = np.asarray(sv["cumulative_times"], dtype=float)

    events = []
    admix_i = 0
    for e, ev in enumerate(dem.ordered_events):
        t_mean = float(cum[:, e].mean())
        t_sd = float(cum[:, e].std())
        if ev["type"] == "MERGE":
            rec = {"event": e + 1, "type": "MERGE", "node": ev["parent"],
                   "children": list(ev["children"]),
                   "time_mean": t_mean, "time_sd": t_sd}
        else:  # ADMIXTURE
            child = ev["child"]
            p1, p2 = ev["parents"]
            f_mean, f_sd = col("admixture_fractions", admix_i)
            rec = {"event": e + 1, "type": "ADMIXTURE", "node": child,
                   "parents": [p1, p2],
                   "time_mean": t_mean, "time_sd": t_sd,
                   "frac_to_parent1": f_mean, "frac_to_parent1_sd": f_sd,
                   "parent1": p1}
            admix_i += 1
        events.append(rec)

    # Effective sizes (haploid).
    ne = {}
    if varying:
        Ne = np.asarray(sv["Ne"], dtype=float)
        for k, nm in enumerate(node_names):
            ne[nm] = {"Ne_haploid_mean": float(Ne[:, k].mean()),
                      "Ne_haploid_sd": float(Ne[:, k].std())}
    else:
        eff = np.asarray(sv["effective_N"], dtype=float)
        shared = {"Ne_haploid_mean": float(eff.mean()),
                  "Ne_haploid_sd": float(eff.std())}
        for nm in node_names:
            ne[nm] = shared

    out = {"events": events, "effective_sizes": ne}
    if "kappa_snp" in sv:
        k = np.asarray(sv["kappa_snp"], dtype=float)
        out["kappa_snp"] = {"mean": float(k.mean()), "sd": float(k.std())}
    return out


def apply_estimates_to_topology(dem, est):
    """Write fitted times / fractions / Ne onto the topology for plotting."""
    for rec in est["events"]:
        if rec["type"] == "MERGE":
            dem.set_merge_time(rec["node"], rec["time_mean"])
        else:
            dem.set_admixture_parameters(
                rec["node"], time=rec["time_mean"],
                fraction_parent_1=rec["frac_to_parent1"],
                parent_1_name=rec["parent1"],
            )
    dem.finalize_root()
    # Store Ne as DIPLOID (codebase convention: haploid = 2 * node.ne).
    for nm, d in est["effective_sizes"].items():
        dem.set_node_ne(nm, d["Ne_haploid_mean"] / 2.0)


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------
def print_report(dem, est, varying):
    print("\n" + "=" * 72)
    print("PARAMETERS ALONG THE TOPOLOGY  (posterior mean +/- sd)")
    print("=" * 72)

    print(f"\n{'Event':>5} {'Type':<10} {'Node':<10} "
          f"{'Time (gen)':>20} {'Admix frac (->parent1)':>26}")
    print("-" * 72)
    for rec in est["events"]:
        t = f"{rec['time_mean']:.1f} +/- {rec['time_sd']:.1f}"
        if rec["type"] == "ADMIXTURE":
            f = (f"{rec['frac_to_parent1']:.3f} +/- "
                 f"{rec['frac_to_parent1_sd']:.3f} (->{rec['parent1']})")
        else:
            f = "-"
        print(f"{rec['event']:>5} {rec['type']:<10} {rec['node']:<10} "
              f"{t:>20} {f:>26}")

    print(f"\n{'Node':<10} {'Ne (haploid)':>24}")
    print("-" * 36)
    if varying:
        for nm, d in est["effective_sizes"].items():
            print(f"{nm:<10} {d['Ne_haploid_mean']:>12.0f} +/- "
                  f"{d['Ne_haploid_sd']:>8.0f}")
    else:
        d = next(iter(est["effective_sizes"].values()))
        print(f"{'(shared)':<10} {d['Ne_haploid_mean']:>12.0f} +/- "
              f"{d['Ne_haploid_sd']:>8.0f}")

    if "kappa_snp" in est:
        print(f"\nkappa_snp: {est['kappa_snp']['mean']:.3f} +/- "
              f"{est['kappa_snp']['sd']:.3f}")
    if "_meta" in est and "cl_power" in est["_meta"]:
        print(f"\ncl_power (fixed global composite-likelihood power): "
              f"{est['_meta']['cl_power']:g}")


# ----------------------------------------------------------------------
# Optional data generation (shell out to generate_stan_data.py)
# ----------------------------------------------------------------------
def generate_data(args):
    """Run generate_stan_data.py to build <args.data>.npz/.json."""
    cmd = [
        sys.executable, os.path.join(_THIS_DIR, "generate_stan_data.py"),
        "--out", args.data,
        "--ibd-method", args.ibd_method,
        "--pop-order", *args.pop_order,
    ]
    if args.chrom:
        cmd += ["--chrom", str(args.chrom)]
    if args.data_dir:
        cmd += ["--data-dir", args.data_dir]
    if args.ibd_glob:
        cmd += ["--ibd-glob", args.ibd_glob]
    print(f"[generate] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Fit a topology to real SNP/IBD data and report the "
                    "parameters along the topology."
    )
    ap.add_argument("--data", default=os.path.join(_THIS_DIR, "real_stan_data"),
                    help="Data product prefix (<prefix>.npz + <prefix>.json).")
    ap.add_argument("--topology", default="amr_admixture",
                    help="Built-in name (amr_admixture, tree_no_admix) or path "
                         "to a .py defining build_topology(pop_order).")
    ap.add_argument("--model", choices=["mixed", "ibd", "snp"], default="mixed")
    ap.add_argument("--nmodel", choices=["Nvarying", "Nfixed"],
                    default="Nvarying", help="Per-node Ne (Nvarying) or shared "
                                             "effective_N (Nfixed).")
    ap.add_argument("--pop-order", nargs="+", default=DEFAULT_POP_ORDER,
                    help="Must match the pop_order saved in the data .json.")
    ap.add_argument("--t-max", type=float, default=DEFAULT_T_MAX,
                    help="Upper bound for the final open IBD epoch (generations).")
    ap.add_argument("--se-floor", type=float, default=1e-8,
                    help="Minimum IBD standard error. The default 1e-8 is tiny "
                         "next to real ibd_hat (~1e-5); raising it (e.g. 1e-6) "
                         "softens the over-sharp IBD Gaussian that can cause "
                         "non-finite gradients on real data.")
    ap.add_argument("--cl-power", type=float, default=1e-3,
                    help="Fixed global composite-likelihood power c in (0,1] "
                         "(mixed/ibd models). Tempers every independent log-density "
                         "term to correct composite-likelihood overconfidence. A "
                         "learned c rails to 0, so it is fixed/profiled here.")
    ap.add_argument("--out", default=None,
                    help="Output prefix (default: <data>_fit).")
    # --- optional end-to-end data generation ---
    ap.add_argument("--generate", action="store_true",
                    help="Build the data product first via generate_stan_data.py.")
    ap.add_argument("--ibd-method", choices=["hapibd", "fastsmc", "file", "none"],
                    default="hapibd",
                    help="IBD caller for --generate (default hapibd).")
    ap.add_argument("--chrom", default=None,
                    help="Restrict generation to one chromosome (e.g. 22).")
    ap.add_argument("--data-dir", default=None,
                    help="merged_pruned dir for --generate (default in script).")
    ap.add_argument("--ibd-glob", default=None,
                    help="With --ibd-method file: glob of precomputed segments.")
    args = ap.parse_args()

    varying = args.nmodel == "Nvarying"
    out_prefix = args.out or (args.data + "_fit")

    # 1. (optional) build the SNP/IBD product
    if args.generate:
        generate_data(args)

    npz_path = args.data + ".npz"
    json_path = args.data + ".json"
    for p in (npz_path, json_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"{p} not found. Build it with generate_stan_data.py "
                f"(or pass --generate)."
            )

    # 2. topology
    dem = load_topology(args.topology, args.pop_order)
    dem.is_valid()
    print(f"\n[topology] {args.topology}: {len(dem.ordered_events)} events, "
          f"{dem.n_admix} admixture(s), {len(dem.nodes)} nodes")

    # 3. assemble Stan data (reuses the build_*_stan_data path)
    stan_data = assemble_stan_data(
        dem, npz_path, json_path, model=args.model, T_max=args.t_max,
    )
    if args.model in ("mixed", "ibd"):
        stan_data["cl_power"] = args.cl_power
        print(f"[cl_power] fixed composite-likelihood power c = {args.cl_power:g}")

    # 4. compile + Pathfinder fit
    stan_path = stan_file_for(args.model, varying)
    print(f"[stan] compiling {os.path.basename(stan_path)}")
    model = CmdStanModel(stan_file=stan_path)
    init = build_init(dem, args.model, varying)
    print("[stan] running Pathfinder ...")
    fit = model.pathfinder(data=stan_data, inits=init, psis_resample=True,
                           show_console=False)
    sv = fit.stan_variables()

    # 5. map estimates onto the topology, report + save
    est = extract_estimates(dem, sv, varying)
    est["_meta"] = {"topology": args.topology, "model": args.model,
                    "nmodel": args.nmodel, "data": args.data,
                    "pop_order": list(args.pop_order), "t_max": args.t_max,
                    "cl_power": args.cl_power}
    print_report(dem, est, varying)

    est_path = out_prefix + "_estimates.json"
    with open(est_path, "w") as f:
        json.dump(est, f, indent=2)
    print(f"\n[save] {est_path}")

    apply_estimates_to_topology(dem, est)
    fig, ax = plt.subplots(figsize=(10, 6))
    dem.plot_demography(scale=True, ax=ax)
    ax.set_title(f"{args.topology} fitted ({args.model}, {args.nmodel})")
    fig_path = out_prefix + "_topology.png"
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fig_path}")


if __name__ == "__main__":
    main()
