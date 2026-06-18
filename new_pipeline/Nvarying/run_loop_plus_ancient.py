"""
Loop-plus-ancient-admixture experiment:
  Constructs a 4-population scenario where exactly one event is invisible
  to each modality:

    - A "loop" on lineage a (recent): a splits backward into aP1 (Ne=Nx) and
      aP2 (Ne=Ny), drifting independently for Δt₂, then merges back into
      aTop. With α=0.5, Nx=Ny=N/2, the per-pair coalescent rate during the
      loop is α²/Nx + (1-α)²/Ny = 1/N — drift contribution matches a single
      population of size N exactly to first order. f-statistics see only
      the integrated drift, so the loop is f-stat-INVISIBLE.
      IBD sees the bimodal TMRCA distribution: same-side pairs (mass
      α²+(1-α)²) coalesce at rate 1/Nx (or 1/Ny), faster than 1/N, while
      cross-side pairs (mass 2α(1-α)) cannot coalesce within Δt₂ at all.
      => IBD distinguishes loop from no-loop.

    - An ancient admixture on c (deep): c splits backward into cP1 (toward
      ab side, fraction f_OLD) and cP2 (toward d side, 1-f_OLD), at time
      t_OLD. Mean tract length is 100/t_OLD ~0.14 cM, well below the
      0.5 cM IBD bin floor. f-statistics scale as f²·(drift between
      sources) and remain detectable.
      => SNP distinguishes admix from no-admix; IBD does not.

  Topology:                                root
                                          /    \
                                       left    right
                                       / \      / \
                                     ab  cP1   d  cP2
                                    /  \   \____c____/
                                  aTop  b
                                  /  \
                                aP1  aP2  (loop, recent)
                                  \  /
                                   a

  T_true: loop + ancient admix
  T_alt1: loop only (c is a normal sister of (ab,d))
  T_alt2: no loop (a is a single Ne leaf), ancient admix kept

  Predictions:
    IBD-only: dELBO(true-alt1) ~ 0   (blind to ancient admix)
              dELBO(true-alt2) >  0  (sees the loop's bimodal TMRCA)
    SNP-only: dELBO(true-alt1) >  0  (sees the ancient admix in f-stats)
              dELBO(true-alt2) ~ 0   (blind to the loop by construction)
    Mixed   : positive on BOTH axes.
"""

import sys, os, time
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import re

from simulation_methods import (
    simulate_msprime,
    build_ibd_stan_data,
    build_snp_stan_data,
    build_mixed_stan_data,
)
from inference_methods import (
    prepare_snp_blocks,
    pool_snp_blocks,
    resample_snp_covariance,
)
from ibd_jackknife import (
    calculate_ibd_blocks_mrca,
    pool_multiple_simulations,
    resample_ibd_with_jackknife_variance,
)
from demography import DemographicTopology
from relative_error import plot_relative_error_boxplot
from cmdstanpy import CmdStanModel


def extract_elbo(fit):
    best = -np.inf
    for stdout_file in fit._runset.stdout_files:
        try:
            with open(stdout_file, 'r') as f:
                for line in f:
                    m = re.search(r'Best Iter:.*?ELBO \(([-+\d.eE]+)\)', line)
                    if m:
                        v = float(m.group(1))
                        if v > best:
                            best = v
        except (OSError, ValueError):
            continue
    return best if np.isfinite(best) else np.nan


def get_admix_event_mapping(dem):
    mapping = []
    admix_counter = 0
    for i, ev in enumerate(dem.ordered_events):
        if ev['type'] == 'ADMIXTURE':
            mapping.append((ev['child'], i, admix_counter))
            admix_counter += 1
    return mapping


def get_merge_event_index(dem, parent_name):
    """Return the cumulative_times index of the merge event that creates `parent_name`."""
    for i, ev in enumerate(dem.ordered_events):
        if ev['type'] == 'MERGE' and ev['parent'] == parent_name:
            return i
    return None


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 100
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15}
# Haploid sample count per population (diploid samples x ploidy 2)
N_HAPLOID_PER_POP = {_p: 2 * _c for _p, _c in SAMPLES_PER_POP.items()}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 200, 300, 500, 750, 1000]

# Bin floor 0.2 cM extends visibility deep into the loop's TMRCA
# window: with L = 100/(2t), bins down to 0.2 cM capture coalescences
# up to t ~ 250, covering the full loop window [T_LOOP_SPLIT,
# T_LOOP_CLOSE]. Trade-off: ancient admix tracts (mean ~ 0.14 cM) now
# leak ~6% of their mass above 0.2 cM, so IBD-only is no longer
# strictly blind to the ancient admix.
bins = [
    # Small bins (0.2 - 0.5 cM)
    [0.2, 0.25], [0.25, 0.3], [0.3, 0.35], [0.35, 0.4],
    [0.4, 0.45], [0.45, 0.5],
    # Standard bins (0.5 - 50 cM)
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7],
    [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],
    [8.0, 20.0], [20.0, BLOCK_SIZE_CM]
]

# --- Population sizes (haploid; demography API takes diploid = haploid//2) ---
# To make the loop detectable, drift inside the loop must accumulate
# substantially within the IBD-visible TMRCA range. This needs a-lineage
# Ne small. For the SNP signal of c's admix, c's leaf branch must drift
# enough to give detectable f-statistics — same trick as v3 (small Ne_d).
NE_BASE   = 10000   # most internal branches
NE_A_LEAF = 2000    # a leaf and aTop (above the loop close)
NE_LOOP   = 1000    # aP1, aP2 individually -> Nx=Ny=N/2 with N=NE_A_LEAF
NE_C_LEAF = 5000    # c leaf branch (drift amplifier for ancient admix)

# --- Loop on a ---
# Loop window kept inside the IBD-visible TMRCA range:
# L = 100/(2t) > 0.2 cM  <=>  t < 250, so window [5, 100] is fully visible.
T_LOOP_SPLIT = 5
T_LOOP_CLOSE = 100           # Δt₂ = 95
ALPHA_LOOP   = 0.5            # symmetric -> max IBD signal

# --- Ancient admix on c ---
T_OLD = 900
F_OLD = 0.70                  # cP1 fraction (toward ab side)

# --- Other event times ---
T_AB_MERGE   = 300
T_LEFT_MERGE = 1100
T_RIGHT_MERGE= 1300
T_ROOT       = 1500

SNP_CUTOFF_TIME = T_ROOT
SNP_MIN_MAF = 0.05


# ====================================================================
# 1. Topology builders
# ====================================================================
def _set_uniform_ne(dem, ne_haploid):
    for n in dem.nodes:
        dem.set_node_ne(n, ne_haploid // 2)


def build_t_true():
    """Loop on a + ancient admix on c."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    # Loop on a
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP1', 'aP2', 'aTop')
    # a-side joins b
    d.add_merge_event('aTop', 'b', 'ab')
    # Ancient admix on c
    d.add_admixture_event('c', 'cP1', 'cP2')
    # cP1 joins ab side; cP2 joins d side
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d', 'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')

    _set_uniform_ne(d, NE_BASE)
    d.set_node_ne('a',    NE_A_LEAF // 2)
    d.set_node_ne('aTop', NE_A_LEAF // 2)
    d.set_node_ne('aP1',  NE_LOOP   // 2)
    d.set_node_ne('aP2',  NE_LOOP   // 2)
    d.set_node_ne('c',    NE_C_LEAF // 2)

    d.set_admixture_parameters('a', time=T_LOOP_SPLIT,
                               fraction_parent_1=ALPHA_LOOP,
                               parent_1_name='aP1')
    d.set_merge_time('aTop',  T_LOOP_CLOSE)
    d.set_merge_time('ab',    T_AB_MERGE)
    d.set_admixture_parameters('c', time=T_OLD,
                               fraction_parent_1=F_OLD,
                               parent_1_name='cP1')
    d.set_merge_time('left',  T_LEFT_MERGE)
    d.set_merge_time('right', T_RIGHT_MERGE)
    d.set_merge_time('root',  T_ROOT)
    d.finalize_root()
    return d


def build_t_alt1():
    """Loop on a, NO ancient admix on c. c is a normal lineage joining the
    deep tree as a sister of (ab, d)."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP1', 'aP2', 'aTop')
    d.add_merge_event('aTop', 'b', 'ab')
    d.add_merge_event('ab', 'c', 'abc')
    d.add_merge_event('abc', 'd', 'root')
    d.finalize_root()
    return d


def build_t_alt2():
    """No loop on a (a is a single-Ne leaf), ancient admix on c kept."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_merge_event('a', 'b', 'ab')
    d.add_admixture_event('c', 'cP1', 'cP2')
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d', 'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    d.finalize_root()
    return d


print("=" * 60)
print("Defining topologies")
print("=" * 60)
dem_true = build_t_true()
dem_alt1 = build_t_alt1()
dem_alt2 = build_t_alt2()

topology_dems = [dem_true, dem_alt1, dem_alt2]
topology_labels = [
    "T_true (loop + ancient)",
    "T_alt1 (loop only)",
    "T_alt2 (ancient only)",
]
for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, "
          f"n_admix={d.n_admix}, n_nodes={len(d.nodes)}")

n_topos = len(topology_dems)
dem_sim = dem_true


# ====================================================================
# 2. Simulate and pool blocks
# ====================================================================
packed_list = []
n_blocks_list = []
snp_blocks_list = []
pair_info = None
pop_ids = None

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    t0 = time.time()
    mts = simulate_msprime(
        dem_sim,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )
    t1 = time.time()
    print(f"  [msprime] {t1 - t0:.2f}s  |  "
          f"n_trees={mts.num_trees}  n_sites={mts.num_sites}")

    packed, n_blocks, pop_samples, pids, pi = calculate_ibd_blocks_mrca(
        mts, bins=bins,
        block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
    )
    t2 = time.time()
    print(f"  [IBD] {t2 - t1:.2f}s")
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi
        pop_ids = pids

    snp_blks, _ = prepare_snp_blocks(
        mts, dem_sim,
        block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME, min_maf=SNP_MIN_MAF,
    )
    t3 = time.time()
    print(f"  [SNP] {t3 - t2:.2f}s")
    snp_blocks_list.append(snp_blks)

pooled_ibd, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
pooled_snp = pool_snp_blocks(snp_blocks_list)

print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} simulations")


# ====================================================================
# 3. Compile Stan models
# ====================================================================
ibd_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nvarying.stan"))
snp_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "snp_model_Nvarying.stan"))
mixed_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nvarying.stan"))

model_labels = ["IBD-only", "SNP-only", "Mixed"]
model_stans = [ibd_stan, snp_stan, mixed_stan]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]


# ====================================================================
# 4. Per-topology event indices for parameter extraction
# ====================================================================
# Indices into cumulative_times for the events whose timing we recover.
# loop_split = a's ADMIXTURE event; loop_close = aTop merge event.
# c_admix = c's ADMIXTURE event.
admix_mappings = [get_admix_event_mapping(d) for d in topology_dems]
loop_close_idx = [get_merge_event_index(d, 'aTop') for d in topology_dems]

for ti, (lbl, mp, lci) in enumerate(zip(topology_labels, admix_mappings, loop_close_idx)):
    print(f"  {lbl} admix mapping: {mp}, loop_close cum_t idx: {lci}")


# ====================================================================
# 5. Run replicates
# ====================================================================
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}
admix_params = {
    label: {ti: {cm: [] for cm in cm_values} for ti in range(n_topos)}
    for label in model_labels
}

# Per-replicate stan_variables() for the relative-error box plot:
# Mixed model, true topology (ti == 0), largest genome length only.
rel_err_stanvars = []
REL_ERR_CM = cm_values[-1]

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} reps x {n_topos} topos x {len(model_labels)} models)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed, replace=False)

        ibd_mean, ibd_var = resample_ibd_with_jackknife_variance(
            pooled_ibd, n_blocks_total=total_blocks,
            pop_ids=pop_ids, pair_info=pair_info, bins=bins,
            target_cm=cm_val, block_size_cm=BLOCK_SIZE_CM,
            chosen_blocks=chosen_blocks,
        )

        leaves = dem_sim.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid, se_block_size=50,
        )

        for m_idx, (m_label, m_stan) in enumerate(zip(model_labels, model_stans)):
            rep_elbos = [np.nan] * n_topos

            for ti, dem_infer in enumerate(topology_dems):
                n_events_t = len(dem_infer.ordered_events)
                n_nodes_t = len(dem_infer.nodes)
                n_admix_t = dem_infer.n_admix

                init = {
                    "times": [100.0] * n_events_t,
                    "mu_log": float(np.log(15000.0)),
                    "sigma_log": 0.3,
                    "Ne_raw": [0.0] * n_nodes_t,
                }
                if n_admix_t > 0:
                    init["admixture_fractions"] = [0.5] * n_admix_t
                if m_label == "Mixed":
                    init["kappa_snp"] = 1.0

                try:
                    if m_label == "IBD-only":
                        sd = build_ibd_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, n_samples_per_pop=N_HAPLOID_PER_POP,
                            T_max=100000, cm=cm_val,
                        )
                    elif m_label == "SNP-only":
                        sd = build_snp_stan_data(
                            dem_infer, w_hat, w_se,
                        )
                    else:
                        sd = build_mixed_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se, n_samples_per_pop=N_HAPLOID_PER_POP,
                            T_max=100000, cm=cm_val,
                        )

                    fit = m_stan.pathfinder(
                        data=sd, inits=init, show_console=False,
                        psis_resample=True,
                    )
                    rep_elbos[ti] = extract_elbo(fit)

                    all_vars = fit.stan_variables()

                    # Collect for the relative-error box plot.
                    if (m_label == "Mixed" and ti == 0
                            and cm_val == REL_ERR_CM):
                        rel_err_stanvars.append(all_vars)

                    cum_t = all_vars['cumulative_times']
                    admix_f = all_vars.get('admixture_fractions', None)
                    pdict = {}
                    for child, ct_idx, af_idx in admix_mappings[ti]:
                        pdict[f'{child}_time'] = cum_t[:, ct_idx].mean()
                        if admix_f is not None:
                            if admix_f.ndim > 1:
                                pdict[f'{child}_frac'] = admix_f[:, af_idx].mean()
                            else:
                                pdict[f'{child}_frac'] = float(admix_f.mean())
                    if loop_close_idx[ti] is not None:
                        pdict['aTop_time'] = cum_t[:, loop_close_idx[ti]].mean()
                    admix_params[m_label][ti][cm_val].append(pdict)

                except Exception as e:
                    if rep < 3:
                        print(f"    [WARN] {m_label} {topology_labels[ti]} "
                              f"rep {rep+1}: {e}")

            if all(np.isfinite(rep_elbos)):
                elbo_results[m_label][cm_val].append(tuple(rep_elbos))


# ====================================================================
# 6a. Plot 1: dELBO comparisons
# ====================================================================
fig, axes = plt.subplots(2, 3, figsize=(28, 16))
x = np.arange(len(cm_values))

comparisons = [
    ("T_true - T_alt1 (detects ANCIENT admix?)", 0, 1),
    ("T_true - T_alt2 (detects LOOP?)",          0, 2),
]

for row_idx, (comp_label, idx_a, idx_b) in enumerate(comparisons):
    for m_idx, m_label in enumerate(model_labels):
        color = model_colors[m_idx]
        ax = axes[row_idx, m_idx]

        diff_data = []
        for cm_val in cm_values:
            elbos = elbo_results[m_label][cm_val]
            if len(elbos) > 0:
                arr = np.array(elbos)
                diff_data.append(arr[:, idx_a] - arr[:, idx_b])
            else:
                diff_data.append(np.array([np.nan]))

        vp = ax.violinplot(
            diff_data, positions=x, widths=0.6,
            showmeans=True, showmedians=True, showextrema=False,
        )
        for body in vp['bodies']:
            body.set_facecolor(color)
            body.set_alpha(0.5)
            body.set_edgecolor("k")
            body.set_linewidth(0.5)
        vp['cmeans'].set_color("k")
        vp['cmeans'].set_linewidth(1)
        vp['cmedians'].set_color("k")
        vp['cmedians'].set_linewidth(1)
        vp['cmedians'].set_linestyle("--")

        ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=10)
        ax.set_xlabel("Genome length (cM)", fontsize=11)
        ax.set_ylabel("delta-ELBO", fontsize=11)
        ax.set_title(f"{m_label}: {comp_label}", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.15)

        for ci, cm_val in enumerate(cm_values):
            dd = diff_data[ci]
            if not np.all(np.isnan(dd)):
                pct_win = 100.0 * np.mean(dd > 0)
                ax.text(ci, ax.get_ylim()[0], f"{pct_win:.0f}%",
                        ha='center', va='bottom', fontsize=8,
                        color='green', fontweight='bold')

fig.suptitle(
    f"Loop + ancient admix: T_true vs T_alt1 (loop only) vs T_alt2 (ancient only)\n"
    f"loop on a: t={T_LOOP_SPLIT}->{T_LOOP_CLOSE}, alpha={ALPHA_LOOP}, "
    f"Nx=Ny={NE_LOOP} (matches single Ne={NE_A_LEAF})  |  "
    f"ancient on c: t={T_OLD}, f={F_OLD}, Ne_c={NE_C_LEAF}  |  "
    f"bins from {bins[0][0]} cM",
    fontsize=12, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("loop_plus_ancient_dELBO.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: loop_plus_ancient_dELBO.png")


# ====================================================================
# 6b. Plot 2: Raw ELBO per topology per model
# ====================================================================
TOPO_COLORS_ALL = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}
TOPO_NAMES_ALL  = {0: "T_true", 1: "T_alt1", 2: "T_alt2"}

fig_e, axes_e = plt.subplots(1, 3, figsize=(28, 8))
n_g = n_topos
bar_width = 0.8 / n_g

for m_idx, m_label in enumerate(model_labels):
    ax = axes_e[m_idx]

    for ti in range(n_topos):
        color = TOPO_COLORS_ALL[ti]
        data_list, positions = [], []
        for ci, cm_val in enumerate(cm_values):
            elbos = elbo_results[m_label][cm_val]
            if len(elbos) > 0:
                arr = np.array(elbos)
                data_list.append(arr[:, ti])
            else:
                data_list.append(np.array([np.nan]))
            positions.append(ci + (ti - (n_g - 1) / 2) * bar_width)

        vp = ax.violinplot(
            data_list, positions=positions,
            widths=bar_width * 0.85,
            showmeans=True, showmedians=False, showextrema=False,
        )
        for body in vp['bodies']:
            body.set_facecolor(color)
            body.set_alpha(0.55)
            body.set_edgecolor("k")
            body.set_linewidth(0.5)
        vp['cmeans'].set_color(color)
        vp['cmeans'].set_linewidth(1.5)

    ax.set_xticks(range(len(cm_values)))
    ax.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=10)
    ax.set_xlabel("Genome length (cM)", fontsize=11)
    ax.set_ylabel("ELBO", fontsize=11)
    ax.set_title(f"{m_label}: raw ELBO by topology", fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.15)

    handles = [Patch(facecolor=TOPO_COLORS_ALL[ti], alpha=0.55,
                     label=TOPO_NAMES_ALL[ti]) for ti in range(n_topos)]
    ax.legend(handles=handles, fontsize=9, loc='best')

fig_e.suptitle("Raw ELBO per topology per model (higher = better fit)",
               fontsize=13, fontweight="bold")
fig_e.tight_layout(rect=[0, 0, 1, 0.94])
fig_e.savefig("loop_plus_ancient_ELBO.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: loop_plus_ancient_ELBO.png")


# ====================================================================
# 6c. Plot 3: Parameter recovery
#     Loop params: compare T_true vs T_alt1 (both have loop)
#     Ancient params: compare T_true vs T_alt2 (both have ancient admix)
# ====================================================================
TOPO_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}

param_rows = [
    ("Loop split time (a admix)",  "a_time",    float(T_LOOP_SPLIT),
     [(0, "T_true"), (1, "T_alt1")]),
    ("Loop close time (aTop merge)", "aTop_time", float(T_LOOP_CLOSE),
     [(0, "T_true"), (1, "T_alt1")]),
    ("Loop fraction (a frac)",     "a_frac",    ALPHA_LOOP,
     [(0, "T_true"), (1, "T_alt1")]),
    ("Ancient admix time (c)",     "c_time",    float(T_OLD),
     [(0, "T_true"), (2, "T_alt2")]),
    ("Ancient admix frac (c)",     "c_frac",    F_OLD,
     [(0, "T_true"), (2, "T_alt2")]),
]

fig2, axes2 = plt.subplots(len(param_rows), 3, figsize=(24, 5 * len(param_rows)))

for row_idx, (row_label, param_key, truth, topo_entries) in enumerate(param_rows):
    for m_idx, m_label in enumerate(model_labels):
        ax = axes2[row_idx, m_idx]
        n_g = len(topo_entries)
        bar_width = 0.7 / n_g

        for gi, (ti, topo_disp) in enumerate(topo_entries):
            color = TOPO_COLORS[ti]
            data_list, positions = [], []
            for ci, cm_val in enumerate(cm_values):
                vals = [p[param_key] for p in admix_params[m_label][ti][cm_val]
                        if param_key in p]
                data_list.append(np.array(vals) if len(vals) > 0
                                 else np.array([np.nan]))
                positions.append(ci + (gi - (n_g - 1) / 2) * bar_width)

            vp = ax.violinplot(
                data_list, positions=positions,
                widths=bar_width * 0.85,
                showmeans=True, showmedians=False, showextrema=False,
            )
            for body in vp['bodies']:
                body.set_facecolor(color)
                body.set_alpha(0.55)
                body.set_edgecolor("k")
                body.set_linewidth(0.5)
            vp['cmeans'].set_color(color)
            vp['cmeans'].set_linewidth(1.5)

        ax.axhline(truth, color='red', ls='--', lw=1.2,
                   label=f"truth = {truth:g}")
        ax.set_xticks(range(len(cm_values)))
        ax.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=9)
        ax.set_xlabel("Genome length (cM)", fontsize=10)
        ax.set_ylabel(param_key, fontsize=10)
        ax.set_title(f"{m_label}: {row_label}", fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.15)

        if m_idx == 0:
            handles = [Patch(facecolor=TOPO_COLORS[ti], alpha=0.55, label=name)
                       for ti, name in topo_entries]
            handles.append(plt.Line2D([0], [0], color='red', ls='--',
                                       lw=1.2, label="truth"))
            ax.legend(handles=handles, fontsize=8, loc='best')

fig2.suptitle(
    f"Loop + ancient parameter recovery  |  "
    f"loop: t_split={T_LOOP_SPLIT}, t_close={T_LOOP_CLOSE}, alpha={ALPHA_LOOP}  |  "
    f"ancient: t={T_OLD}, f={F_OLD}",
    fontsize=13, fontweight="bold"
)
fig2.tight_layout(rect=[0, 0, 1, 0.96])
fig2.savefig("loop_plus_ancient_params.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: loop_plus_ancient_params.png")


# ====================================================================
# Parameter relative error (glike-style box plot)
# Mixed model, true topology, largest genome length.
# ====================================================================
plot_relative_error_boxplot(
    dem_true, rel_err_stanvars, varying=True,
    outpath="loop_plus_ancient_relative_error.png",
    title=f"Parameter relative error (Mixed, true topology, {REL_ERR_CM} cM)",
)


# ====================================================================
# 7. Summary table
# ====================================================================
print("\n" + "=" * 110)
print("LOOP + ANCIENT - SUMMARY")
print("=" * 110)
print(f"{'cm':>6s} | {'Model':>10s} | "
      f"{'dELBO(true-alt1)':>18s} | {'%win_alt1':>10s} | "
      f"{'dELBO(true-alt2)':>18s} | {'%win_alt2':>10s} | "
      f"{'n_ok':>5s}")
print("-" * 110)

for cm_val in cm_values:
    for m_label in model_labels:
        elbos = elbo_results[m_label][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)
            d1 = arr[:, 0] - arr[:, 1]
            d2 = arr[:, 0] - arr[:, 2]
            pct1 = 100.0 * np.mean(d1 > 0)
            pct2 = 100.0 * np.mean(d2 > 0)
            print(f"{cm_val:6d} | {m_label:>10s} | "
                  f"{np.mean(d1):+17.1f} | {pct1:9.1f}% | "
                  f"{np.mean(d2):+17.1f} | {pct2:9.1f}% | "
                  f"{len(elbos):5d}")
        else:
            print(f"{cm_val:6d} | {m_label:>10s} | "
                  f"{'N/A':>18s} | {'N/A':>10s} | "
                  f"{'N/A':>18s} | {'N/A':>10s} | {0:5d}")
