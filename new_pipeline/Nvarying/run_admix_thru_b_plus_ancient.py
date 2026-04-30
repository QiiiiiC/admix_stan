"""
Admix-through-b plus ancient-admixture experiment:
  Recent admixture on lineage a routed through b (instead of merging
  aP1, aP2 back into a single aTop). The two parents of a's admixture
  follow distinct paths and re-converge only at the deeper ab node:

    - At t0=T_LOOP_SPLIT: a admixes backward into aP1 (frac alpha) and
      aP2 (frac 1-alpha).
    - At t1=T_APB_MERGE: aP2 + b -> bP. So aP2's lineage merges with
      b before reaching ab.
    - At t2=T_LOOP_CLOSE: bP + aP1 -> ab.

    The structure differs from the simple "Nx=Ny=N/2 loop" in two
    important ways:
      1. The aP2 branch passes through b's lineage (via bP). This
         leaks a small amount of f-stat signal into Cov(Y_a, Y_b).
      2. b is observed, so its drift contributions tie into the loop
         parameters -- this REDUCES the parameter degeneracy (you
         can't make Ne_aP1 -> infinity without breaking Var(Y_b)).

    F2(a,b) distortion vs the no-recent-admix alternative:
        F2_diff = alpha^2 * (Δt3 / Ne_bP)
    where Δt3 = t2 - t1 is the bP edge length. With alpha=0.5,
    Δt3=5, Ne_bP=10000, this is ~1.25e-4 -- below typical f-stat
    SE for moderate-genome data, so SNP-only model can't detect it
    cleanly. IBD bimodality remains strong (same-side mass ~5%
    with Ne_aP1=Ne_aP2=1000).

    - An ancient admixture on c (deep): c splits backward into cP1
      (toward ab side, fraction f_OLD) and cP2 (toward d side,
      1-f_OLD), at time t_OLD. Mean tract length is 100/t_OLD ~0.14
      cM. f-statistics scale as f^2 * (drift between sources) and
      remain detectable.
      => SNP distinguishes ancient admix from no-admix; IBD does not.

  Topology:                                root
                                          /    \
                                       left    right
                                       / \      / \
                                     ab  cP1   d  cP2
                                    /  \   \____c____/
                                  aP1   bP
                                  /     /  \
                                  \   aP2   b
                                   \  /
                                    a    (admix at t0; aP2+b->bP at t1)

  T_true: recent admix-through-b + ancient admix on c
  T_alt1: no recent admix (a+b->ab directly), ancient admix kept
  T_alt2: recent admix kept, no ancient admix (c is sister of (ab,d))

  Predictions:
    IBD-only: dELBO(true-alt1) >  0  (sees the bimodal TMRCA)
              dELBO(true-alt2) ~  0  (blind to ancient admix)
    SNP-only: dELBO(true-alt1) ~  0  (recent admix near-invisible by tuning)
              dELBO(true-alt2) >  0  (sees the ancient admix in f-stats)
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

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 300, 500, 750, 1000]

bins = [
    [0.2, 0.25], [0.25, 0.3], [0.3, 0.35], [0.35, 0.4],
    [0.4, 0.45], [0.45, 0.5],
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7],
    [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],
    [8.0, 20.0], [20.0, BLOCK_SIZE_CM]
]

# --- Population sizes (haploid; demography API takes diploid = haploid//2) ---
NE_BASE   = 10000   # most internal branches and b-leaf
NE_A_LEAF = 10000    
NE_LOOP   = 5000    
NE_BP     = 10000   # bP edge (large to minimize F2 distortion)
NE_C_LEAF = 10000    # c leaf branch (drift amplifier for ancient admix)

# --- Recent admix on a (through b) ---
T_LOOP_SPLIT  = 5     # t0: a admixes into aP1, aP2
T_APB_MERGE   = 95    # t1: aP2 + b -> bP
T_LOOP_CLOSE  = 100   # t2: bP + aP1 -> ab          ⇒  Δt3 = 5
ALPHA_LOOP    = 0.5

# --- Ancient admix on c ---
T_OLD = 900
F_OLD = 0.70

# --- Other event times ---
T_LEFT_MERGE  = 1100
T_RIGHT_MERGE = 1300
T_ROOT        = 1500

SNP_CUTOFF_TIME = T_ROOT
SNP_MIN_MAF = 0.05
FIXED_NE = NE_BASE


# ====================================================================
# 1. Topology builders
# ====================================================================
def _set_uniform_ne(dem, ne_haploid):
    for n in dem.nodes:
        dem.set_node_ne(n, ne_haploid // 2)


def build_t_true():
    """Recent admix-through-b on a + ancient admix on c."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    # Recent admix on a, with aP2 routed through b
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP2', 'b',  'bP')      # t1
    d.add_merge_event('bP',  'aP1', 'ab')     # t2
    # Ancient admix on c
    d.add_admixture_event('c', 'cP1', 'cP2')
    # Deep tree
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d',  'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')

    _set_uniform_ne(d, NE_BASE)
    d.set_node_ne('a',    NE_A_LEAF // 2)
    d.set_node_ne('aP1',  NE_LOOP   // 2)
    d.set_node_ne('aP2',  NE_LOOP   // 2)
    d.set_node_ne('bP',   NE_BP     // 2)
    d.set_node_ne('c',    NE_C_LEAF // 2)

    d.set_admixture_parameters('a', time=T_LOOP_SPLIT,
                               fraction_parent_1=ALPHA_LOOP,
                               parent_1_name='aP1')
    d.set_merge_time('bP',    T_APB_MERGE)
    d.set_merge_time('ab',    T_LOOP_CLOSE)
    d.set_admixture_parameters('c', time=T_OLD,
                               fraction_parent_1=F_OLD,
                               parent_1_name='cP1')
    d.set_merge_time('left',  T_LEFT_MERGE)
    d.set_merge_time('right', T_RIGHT_MERGE)
    d.set_merge_time('root',  T_ROOT)
    d.finalize_root()
    return d


def build_t_alt1():
    """No recent admix on a (a + b -> ab directly), ancient admix on c kept."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_merge_event('a', 'b', 'ab')
    d.add_admixture_event('c', 'cP1', 'cP2')
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d',  'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    d.finalize_root()
    return d


def build_t_alt2():
    """Recent admix-through-b kept, NO ancient admix on c (c is a normal
    lineage joining the deep tree as a sister of (ab, d))."""
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP2', 'b',  'bP')
    d.add_merge_event('bP',  'aP1', 'ab')
    d.add_merge_event('ab', 'd',   'abd')
    d.add_merge_event('abd', 'c',  'root')
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
    "T_true (recent thru b + ancient)",
    "T_alt1 (ancient only)",
    "T_alt2 (recent only)",
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
# Track the recent admix event (a) and the two intermediate merges
# (bP = aP2+b, ab = bP+aP1). Ancient admix event on c also tracked.
admix_mappings = [get_admix_event_mapping(d) for d in topology_dems]
bP_close_idx   = [get_merge_event_index(d, 'bP') for d in topology_dems]
ab_close_idx   = [get_merge_event_index(d, 'ab') for d in topology_dems]

for ti, (lbl, mp, bpi, abi) in enumerate(zip(
        topology_labels, admix_mappings, bP_close_idx, ab_close_idx)):
    print(f"  {lbl} admix mapping: {mp}, bP cum_t idx: {bpi}, ab cum_t idx: {abi}")


# ====================================================================
# 5. Run replicates
# ====================================================================
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}
admix_params = {
    label: {ti: {cm: [] for cm in cm_values} for ti in range(n_topos)}
    for label in model_labels
}

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
                    "Ne": [10000.0] * n_nodes_t,
                }
                if n_admix_t > 0:
                    init["admixture_fractions"] = [0.5] * n_admix_t
                if m_label == "Mixed":
                    init["kappa_snp"] = 1.0

                try:
                    if m_label == "IBD-only":
                        sd = build_ibd_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins,
                            T_max=100000, cm=cm_val,
                        )
                    elif m_label == "SNP-only":
                        sd = build_snp_stan_data(
                            dem_infer, w_hat, w_se, effective_N=FIXED_NE,
                        )
                    else:
                        sd = build_mixed_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se,
                            T_max=100000, cm=cm_val,
                        )

                    fit = m_stan.pathfinder(
                        data=sd, inits=init, show_console=False,
                        psis_resample=True,
                    )
                    rep_elbos[ti] = extract_elbo(fit)

                    all_vars = fit.stan_variables()
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
                    if bP_close_idx[ti] is not None:
                        pdict['bP_time'] = cum_t[:, bP_close_idx[ti]].mean()
                    if ab_close_idx[ti] is not None:
                        pdict['ab_time'] = cum_t[:, ab_close_idx[ti]].mean()
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
    ("T_true - T_alt1 (detects RECENT admix?)", 0, 1),
    ("T_true - T_alt2 (detects ANCIENT admix?)", 0, 2),
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
    f"Recent admix-thru-b + ancient: T_true vs T_alt1 (no recent) vs T_alt2 (no ancient)\n"
    f"recent on a: t={T_LOOP_SPLIT}->{T_APB_MERGE}->{T_LOOP_CLOSE}, alpha={ALPHA_LOOP}, "
    f"Ne_aP1=Ne_aP2={NE_LOOP}, Ne_bP={NE_BP}  |  "
    f"ancient on c: t={T_OLD}, f={F_OLD}, Ne_c={NE_C_LEAF}  |  "
    f"bins from {bins[0][0]} cM",
    fontsize=12, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("admix_thru_b_plus_ancient_dELBO.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_dELBO.png")


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
fig_e.savefig("admix_thru_b_plus_ancient_ELBO.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_ELBO.png")


# ====================================================================
# 6c. Plot 3: Parameter recovery
#     Recent params: compare T_true vs T_alt2 (both have recent admix)
#     Ancient params: compare T_true vs T_alt1 (both have ancient admix)
# ====================================================================
TOPO_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}

param_rows = [
    ("Recent admix split time (a)", "a_time",  float(T_LOOP_SPLIT),
     [(0, "T_true"), (2, "T_alt2")]),
    ("aP2 + b merge time (bP)",     "bP_time", float(T_APB_MERGE),
     [(0, "T_true"), (2, "T_alt2")]),
    ("bP + aP1 merge time (ab)",    "ab_time", float(T_LOOP_CLOSE),
     [(0, "T_true"), (2, "T_alt2")]),
    ("Recent admix fraction (a)",   "a_frac",  ALPHA_LOOP,
     [(0, "T_true"), (2, "T_alt2")]),
    ("Ancient admix time (c)",      "c_time",  float(T_OLD),
     [(0, "T_true"), (1, "T_alt1")]),
    ("Ancient admix frac (c)",      "c_frac",  F_OLD,
     [(0, "T_true"), (1, "T_alt1")]),
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
    f"Recent thru-b + ancient parameter recovery  |  "
    f"recent: t_split={T_LOOP_SPLIT}, t_apb={T_APB_MERGE}, t_ab={T_LOOP_CLOSE}, alpha={ALPHA_LOOP}  |  "
    f"ancient: t={T_OLD}, f={F_OLD}",
    fontsize=13, fontweight="bold"
)
fig2.tight_layout(rect=[0, 0, 1, 0.96])
fig2.savefig("admix_thru_b_plus_ancient_params.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_params.png")


# ====================================================================
# 7. Summary table
# ====================================================================
print("\n" + "=" * 110)
print("RECENT-THRU-B + ANCIENT - SUMMARY")
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
