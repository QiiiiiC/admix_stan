"""
Relative-error box plots (glike-style supplementary figure).

For a single inference setting (one model, one data/genome-length setting, the
true data-generating topology) we collect, across replicates, the posterior-mean
point estimate of every model parameter and compare it to the known simulation
truth.  The figure is a box plot with

    y-axis : relative error  (estimate - truth) / truth
    x-axis : parameter name

Parameters are, in order:
    * event times      (one per topology event, chronological)
    * admixture fracs  (one per admixture event, chronological)
    * effective sizes  (Nfixed: a single box "Ne";
                        Nvarying: one box per population node "Ne[<node>]")

The figure width scales with the number of parameters so that the Nvarying case
(one Ne per node) stays legible.

Conventions (must match simulation_methods / the Stan models):
    * Stan node index a (1-based) == position of the node in dem.nodes (0-based).
    * Stan effective sizes are HAPLOID; dem stores DIPLOID node.ne
      (set via set_node_ne(name, ne_haploid // 2)) so true haploid = 2 * node.ne.
    * Stan admixture_fractions[i] is the fraction to the FIRST listed parent of
      the i-th admixture event (chronological).
"""

import numpy as np
import matplotlib.pyplot as plt


def _col_mean(arr, idx):
    """Posterior mean of column `idx` of a stan_variables() entry.

    Handles both 2-D (draws, dim) arrays and 1-D (draws,) arrays returned for
    length-1 vector/array parameters.
    """
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        return float(arr.mean()) if idx == 0 else float("nan")
    return float(arr[:, idx].mean())


def build_param_spec(dem_true, varying):
    """Return (names, truths, getters) for every parameter of `dem_true`.

    getters[k](stan_vars) -> the posterior-mean point estimate of parameter k
    from one replicate's fit.stan_variables() dict.
    """
    names, truths, getters = [], [], []
    node_names = list(dem_true.nodes.keys())

    # ---- Event times (chronological; one per ordered event) ----
    for e, ev in enumerate(dem_true.ordered_events):
        if ev["type"] == "MERGE":
            p = ev["parent"]
            true_t = dem_true.nodes[p].time_start
            label = f"t[{p}]"
        else:  # ADMIXTURE
            c = ev["child"]
            true_t = dem_true.nodes[c].time_end
            label = f"t[{c}]*"
        if true_t is None:
            continue
        names.append(label)
        truths.append(float(true_t))
        getters.append(lambda sv, e=e: _col_mean(sv["cumulative_times"], e))

    # ---- Admixture fractions (chronological) ----
    ai = 0
    for ev in dem_true.ordered_events:
        if ev["type"] != "ADMIXTURE":
            continue
        child = ev["child"]
        p1 = ev["parents"][0]
        true_f = dem_true.nodes[child].admixture_fractions.get(p1)
        if true_f is not None:
            names.append(f"f[{child}]")
            truths.append(float(true_f))
            getters.append(
                lambda sv, i=ai: _col_mean(sv["admixture_fractions"], i)
            )
        ai += 1

    # ---- Effective population sizes (haploid) ----
    if varying:
        for k, nm in enumerate(node_names):
            ne = dem_true.nodes[nm].ne
            if ne is None:
                continue
            names.append(f"Ne[{nm}]")
            truths.append(2.0 * float(ne))
            getters.append(lambda sv, k=k: _col_mean(sv["Ne"], k))
    else:
        node_ne = [dem_true.nodes[nm].ne for nm in node_names
                   if dem_true.nodes[nm].ne is not None]
        true_ne = 2.0 * float(np.median(node_ne)) if node_ne else float("nan")
        names.append("Ne")
        truths.append(true_ne)
        getters.append(lambda sv: float(np.mean(np.asarray(sv["effective_N"]))))

    return names, np.asarray(truths, dtype=float), getters


def plot_relative_error_boxplot(dem_true, stan_vars_list, varying,
                                outpath, title=""):
    """Draw and save the relative-error box plot.

    Parameters
    ----------
    dem_true : DemographicTopology
        The data-generating (true) topology.
    stan_vars_list : list of dict
        One fit.stan_variables() dict per replicate, all from the SAME
        (model, genome length, true topology) setting.
    varying : bool
        True for the Nvarying models (per-node Ne), False for Nfixed.
    outpath : str
        Where to save the PNG.
    title : str
        Figure title.
    """
    if len(stan_vars_list) == 0:
        print(f"[relative_error] no replicates collected; skipping {outpath}")
        return None

    names, truths, getters = build_param_spec(dem_true, varying)
    n_par = len(names)
    n_rep = len(stan_vars_list)

    est = np.full((n_rep, n_par), np.nan)
    for r, sv in enumerate(stan_vars_list):
        for k, get in enumerate(getters):
            try:
                est[r, k] = get(sv)
            except Exception:
                est[r, k] = np.nan

    with np.errstate(divide="ignore", invalid="ignore"):
        rel = (est - truths[None, :]) / truths[None, :]

    box_data = [rel[np.isfinite(rel[:, k]), k] for k in range(n_par)]

    # Width scales with the number of parameters (Nvarying has many more boxes).
    width = max(8.0, 0.55 * n_par + 2.0)
    fig, ax = plt.subplots(figsize=(width, 6))

    ax.boxplot(box_data, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#9ecae1", edgecolor="k"),
               medianprops=dict(color="k"))
    ax.axhline(0.0, color="red", ls="--", lw=1.2)
    ax.set_xticks(range(1, n_par + 1))
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Relative error  (estimate - truth) / truth", fontsize=11)
    ax.set_xlabel("Parameter", fontsize=11)
    ax.set_title(title or "Parameter relative error", fontsize=12,
                 fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[relative_error] saved {outpath}  "
          f"({n_par} params, {n_rep} replicates)")
    return outpath


# ---------------------------------------------------------------------------
# Self-test of the truth/extraction logic (no plotting, no Stan needed).
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from demography import DemographicTopology

    d = DemographicTopology(["a", "b", "c", "d"])
    d.add_admixture_event("a", "aP1", "aP2")
    d.add_merge_event("aP2", "b", "bP")
    d.add_merge_event("bP", "aP1", "ab")
    d.add_admixture_event("c", "cP1", "cP2")
    d.add_merge_event("ab", "cP1", "left")
    d.add_merge_event("d", "cP2", "right")
    d.add_merge_event("left", "right", "root")
    for nm in d.nodes:
        d.set_node_ne(nm, 5000)
    d.set_admixture_parameters("a", time=10, fraction_parent_1=0.5,
                               parent_1_name="aP1")
    d.set_merge_time("bP", 95)
    d.set_merge_time("ab", 125)
    d.set_admixture_parameters("c", time=900, fraction_parent_1=0.7,
                               parent_1_name="cP1")
    d.set_merge_time("left", 1100)
    d.set_merge_time("right", 1300)
    d.set_merge_time("root", 1500)
    d.finalize_root()

    names, truths, getters = build_param_spec(d, varying=True)
    print("Nvarying params:", len(names))
    for nm, tv in zip(names, truths):
        print(f"  {nm:12s} truth={tv}")

    names_f, truths_f, _ = build_param_spec(d, varying=False)
    print("Nfixed params:", len(names_f), "->", names_f[-1], truths_f[-1])

    # fake stan_vars: exact-truth draws -> relative error should be ~0
    n_events = len(d.ordered_events)
    true_cum = []
    for ev in d.ordered_events:
        if ev["type"] == "MERGE":
            true_cum.append(d.nodes[ev["parent"]].time_start)
        else:
            true_cum.append(d.nodes[ev["child"]].time_end)
    sv = {
        "cumulative_times": np.tile(np.array(true_cum, float), (50, 1)),
        "admixture_fractions": np.tile(np.array([0.5, 0.7]), (50, 1)),
        "Ne": np.tile(np.array([2 * 5000.0] * len(d.nodes)), (50, 1)),
        "effective_N": np.full(50, 2 * 5000.0),
    }
    _, tr, gs = build_param_spec(d, varying=True)
    est = np.array([g(sv) for g in gs])
    rel = (est - tr) / tr
    print("max |relative error| with exact-truth draws:", np.nanmax(np.abs(rel)))
