"""
Mixed-beats-both experiment v4:
  Exploits genuine SNP-covariance identifiability limits (Maier et al. 2023):
  f-statistics are blind to the NUMBER OF PULSES between a pair of lineages
  when the net ancestry is preserved. IBD segment length distributions
  encode pulse timing, so IBD can distinguish them.

True topology (T_true): 5 populations (a, b, c, d, e) with THREE admix events:
  RECENT (on b): TWO pulses from a-side, summing to net 50% a-ancestry:
      t=5  (merge-back t=6 ): f_1 = 0.25       -> max IBD segment ~ 8.3 cM
      t=80 (merge-back t=81): f_2 = 0.25/0.75  -> max IBD segment ~ 0.62 cM
    Net: 0.25 + 0.75 * 0.333 = 0.50 a-ancestry.
  OLD (on d): t=700, 30% from c-side.

  Backward events:
    t=5:    b   ADMIX -> bP1_r(a-side, f=0.25),  bM(0.75)
    t=6:    a + bP1_r MERGE -> ab           (a-side pulse-back: 1-gen gap)
    t=80:   bM  ADMIX -> bP2_r(a-side, f=0.333), bMc(c-side, 0.667)
    t=81:   ab + bP2_r MERGE -> abP         (a-side pulse-back: 1-gen gap)
    t=82:   c + bMc  MERGE -> bC            (c-side pulse-back: 2-gen gap)
    t=700:  d   ADMIX -> dP1(c-side, f=0.30), dP2(e-side, 0.70)
    t=900:  bC + dP1 MERGE -> cbd
    t=1100: abP + cbd MERGE -> left
    t=1200: e + dP2  MERGE -> right
    t=1500: left + right MERGE -> root

  All admix-source-fragment populations (bP1_r, bP2_r, bMc) merge back
  into their source (a or c) within 1-2 generations -> sources are a and
  c themselves, not phantom sister populations. bC is the c-side common
  ancestor (c together with b's c-side ancestry) and participates as
  a single lineage in the d-clade.

Alternative topologies:
  T_alt1 (1-pulse recent + old): ONE pulse on b at t=40 with f=0.50.
    Net a-ancestry in b is identical to T_true -> SNP covariance predictions
    are (nearly) identical. But T_alt1's segments cap at 50/41 ~ 1.22 cM, so
    bins [1.5, 8] cM are a structural-zero region only T_true can populate.
    => Tests "can the model resolve pulse structure from segment lengths?"
  T_alt2 (2-pulse recent + NO old admix): d is a regular sister of abc.
    Covariance predictions differ at the old-admix axis. IBD sees the same
    (old admix segments ~0.07 cM are below 0.5 cM bin floor).
    => Tests "can the model detect the old admix event?"

Expected outcome:
  IBD-only: dELBO(true-alt1) >> 0 (sees pulse structure via long-segment tail)
            dELBO(true-alt2) ~ 0 (blind to old admix)
  SNP-only: dELBO(true-alt1) ~ 0 (blind to pulse count)
            dELBO(true-alt2) > 0  (sees old admix)
  Mixed   : dELBO(true-alt1) > 0 AND dELBO(true-alt2) > 0
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
    best_elbo = -np.inf
    for stdout_file in fit._runset.stdout_files:
        try:
            with open(stdout_file, 'r') as f:
                for line in f:
                    m = re.search(r'Best Iter:.*?ELBO \(([-+\d.eE]+)\)', line)
                    if m:
                        v = float(m.group(1))
                        if v > best_elbo:
                            best_elbo = v
        except (OSError, ValueError):
            continue
    return best_elbo if np.isfinite(best_elbo) else np.nan


def get_admix_event_mapping(dem):
    mapping = []
    admix_counter = 0
    for i, ev in enumerate(dem.ordered_events):
        if ev['type'] == 'ADMIXTURE':
            mapping.append((ev['child'], i, admix_counter))
            admix_counter += 1
    return mapping


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 100
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15, 'e': 15}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

# IBD bins from 0.5 cM (kills old-admix IBD signal). Long range is finer:
# T_alt1's segment-length cap is at 50/41 ~ 1.22 cM, so bins 1.5-8 cM are
# a structural-zero region for T_alt1 and a positive-signal region for T_true.
bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7],
    [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.2], [1.2, 1.5],
    [1.5, 1.7], [1.7, 2.0],
    [2.0, 3.0], [3.0, 5.0],
    [5.0, 8.0], [8.0, 20.0], [20.0, BLOCK_SIZE_CM]
]

# Ne (haploid)
NE_REST = 10000
NE_B = 10000
NE_D = 5000

# --- Recent admix: two pulses on b summing to net 50% from a-side ---
# Pulse 1 is very recent (t=5, merge t=6) -> long IBD segments up to
# 50/6 ~ 8.3 cM. Pulse 2 is well-separated (t=80) so its segment cloud
# is concentrated in [0.5, 0.6] cM and doesn't overlap pulse 1's tail.
# T_alt1's single pulse at t=40 caps segments at 50/41 ~ 1.22 cM, so
# bins [1.5, 8] cM are a structural-zero region for T_alt1.
T_PULSE_1 = 10
T_PULSE_2 = 80
NET_A_IN_B = 0.50
F_PULSE_1 = 0.25
F_PULSE_2 = (NET_A_IN_B - F_PULSE_1) / (1 - F_PULSE_1)   # ~ 0.333
# Single-pulse equivalent (for T_alt1) matches net ancestry:
T_SINGLE = 40
F_SINGLE = NET_A_IN_B

# Gap between an admix event and the immediate pulse-back merge. We don't
# use 1-gen gaps because Stan's `times` vector has `lower=1` -> 1-gen gaps
# sit right at the boundary and pathfinder has trouble there. 5 gens is
# still effectively a pulse (bMc drifts negligibly with Ne=10000) but is
# numerically well-conditioned.
PULSE_BACK_GAP = 5

# --- Old admix on d ---
T_OLD = 700
F_OLD = 0.30

SNP_CUTOFF_TIME = 1500
SNP_MIN_MAF = 0.05
FIXED_NE = NE_REST


# ====================================================================
# 1. Topology builders
# ====================================================================
def _set_all_ne(dem, default):
    for name in dem.nodes:
        dem.set_node_ne(name, default // 2)


def build_t_true():
    """2-pulse on b + old admix on d.

    Going backward, b receives two a-side pulses (at t=T_PULSE_1 and
    t=T_PULSE_2). After the second pulse, the brief c-side fragment (bMc)
    immediately merges with c at t=T_PULSE_2 + 2 to form bC, the c-side
    common ancestor that participates in the d-clade.
    """
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    # Recent: 2 pulses on b
    d.add_admixture_event('b', 'bP1_r', 'bM')      # pulse 1 at T_PULSE_1
    d.add_merge_event('a', 'bP1_r', 'ab')          # a-side pulse-back: T_PULSE_1 + 1
    d.add_admixture_event('bM', 'bP2_r', 'bMc')    # pulse 2 at T_PULSE_2 (bMc = brief c-side intermediate)
    d.add_merge_event('ab', 'bP2_r', 'abP')        # a-side pulse-back: T_PULSE_2 + 1
    d.add_merge_event('c', 'bMc', 'bC')            # c-side pulse-back: T_PULSE_2 + 2 (bC = c-side MRCA)
    # Old admix on d + the rest of the tree
    d.add_admixture_event('d', 'dP1', 'dP2')       # t=700
    d.add_merge_event('bC', 'dP1', 'cbd')          # t=900
    d.add_merge_event('abP', 'cbd', 'left')        # t=1100
    d.add_merge_event('e', 'dP2', 'right')         # t=1200
    d.add_merge_event('left', 'right', 'root')     # t=1500

    _set_all_ne(d, NE_REST)
    d.set_node_ne('b', NE_B // 2)
    d.set_node_ne('bM', NE_B // 2)
    d.set_node_ne('bMc', NE_B // 2)
    d.set_node_ne('d', NE_D // 2)

    d.set_admixture_parameters('b', time=T_PULSE_1,
                               fraction_parent_1=F_PULSE_1,
                               parent_1_name='bP1_r')
    d.set_merge_time('ab', T_PULSE_1 + PULSE_BACK_GAP)             # a-side pulse-back
    d.set_admixture_parameters('bM', time=T_PULSE_2,
                               fraction_parent_1=F_PULSE_2,
                               parent_1_name='bP2_r')
    d.set_merge_time('abP', T_PULSE_2 + PULSE_BACK_GAP)            # a-side pulse-back
    d.set_merge_time('bC',  T_PULSE_2 + 2 * PULSE_BACK_GAP)        # c-side pulse-back
    d.set_admixture_parameters('d', time=T_OLD,
                               fraction_parent_1=F_OLD,
                               parent_1_name='dP1')
    d.set_merge_time('cbd', 900)
    d.set_merge_time('left', 1100)
    d.set_merge_time('right', 1200)
    d.set_merge_time('root', 1500)
    d.finalize_root()
    return d


def build_t_alt1():
    """Single-pulse on b + old admix on d (inference topology, structure only).

    Mirrors v3's T_true: b's a-side ancestry via one pulse (bP1 returns
    to a), b's c-side ancestry (bC) merges with c.
    """
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    d.add_admixture_event('b', 'bP1', 'bC')
    d.add_merge_event('a', 'bP1', 'ab')
    d.add_merge_event('c', 'bC', 'cb')
    d.add_admixture_event('d', 'dP1', 'dP2')
    d.add_merge_event('cb', 'dP1', 'cbd')
    d.add_merge_event('ab', 'cbd', 'left')
    d.add_merge_event('e', 'dP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    d.finalize_root()
    return d


def build_t_alt2():
    """2-pulse on b, no old admix on d (inference topology, structure only)."""
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    d.add_admixture_event('b', 'bP1_r', 'bM')
    d.add_merge_event('a', 'bP1_r', 'ab')
    d.add_admixture_event('bM', 'bP2_r', 'bC')
    d.add_merge_event('ab', 'bP2_r', 'abP')
    d.add_merge_event('c', 'bC', 'cb')
    d.add_merge_event('abP', 'cb', 'abc')
    d.add_merge_event('abc', 'd', 'abcd')
    d.add_merge_event('abcd', 'e', 'root')
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
    "T_true (2-pulse + old)",
    "T_alt1 (1-pulse + old)",
    "T_alt2 (2-pulse, no old)",
]
for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, "
          f"n_admix={d.n_admix}, n_nodes={len(d.nodes)}")

n_topos = len(topology_dems)

dem_sim = dem_true


# ====================================================================
# 2. Simulate and pool blocks (using true topology)
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
          f"n_trees={mts.num_trees}  n_sites={mts.num_sites}  n_nodes={mts.num_nodes}")

    packed, n_blocks, pop_samples, pids, pi = calculate_ibd_blocks_mrca(
        mts, bins=bins,
        block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
    )
    t2 = time.time()
    print(f"  [IBD calc] {t2 - t1:.2f}s")
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
    print(f"  [SNP prep] {t3 - t2:.2f}s")
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
# 4. Run replicates
# ====================================================================
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}
admix_params = {
    label: {ti: {cm: [] for cm in cm_values} for ti in range(n_topos)}
    for label in model_labels
}
admix_mappings = [get_admix_event_mapping(d) for d in topology_dems]
for ti, (lbl, mp) in enumerate(zip(topology_labels, admix_mappings)):
    print(f"  {lbl} admix mapping: {mp}")

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
                    admix_params[m_label][ti][cm_val].append(pdict)

                except Exception as e:
                    if rep < 3:
                        print(f"    [WARN] {m_label} {topology_labels[ti]} "
                              f"rep {rep+1}: {e}")

            if all(np.isfinite(rep_elbos)):
                elbo_results[m_label][cm_val].append(tuple(rep_elbos))


# ====================================================================
# 5a. Plot: ELBO comparisons (2x3 grid)
# ====================================================================
fig, axes = plt.subplots(2, 3, figsize=(28, 16))
x = np.arange(len(cm_values))

comparisons = [
    ("T_true - T_alt1 (detects PULSE STRUCTURE?)", 0, 1),
    ("T_true - T_alt2 (detects OLD admix?)",       0, 2),
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
    f"Mixed beats both v4: 2-pulse recent admix (f1={F_PULSE_1} at t={T_PULSE_1}, "
    f"f2={F_PULSE_2:.2f} at t={T_PULSE_2}, net={NET_A_IN_B}) + "
    f"old admix (t={T_OLD}, f={F_OLD})  |  "
    f"Ne: b={NE_B}, d={NE_D}, rest={NE_REST}  |  bins from {bins[0][0]} cM\n"
    "T_alt1 = 1 pulse (net ancestry matched, SNP-equivalent); "
    "T_alt2 = no old admix (IBD-equivalent given 0.5 cM floor)",
    fontsize=12, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("mixed_beats_both_v4.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: mixed_beats_both_v4.png")


# ====================================================================
# 5b. Plot: Old-admix parameter estimation (2x3 grid)
#     Only the OLD admix is directly comparable across topologies
#     (T_true and T_alt1 both have d admix; T_alt2 does not).
# ====================================================================
TOPO_COLORS = {0: "#1f77b4", 1: "#ff7f0e"}
TOPO_NAMES = {0: "T_true", 1: "T_alt1"}

param_rows = [
    ("Old admix time (d)", "d_time", float(T_OLD), [(0, "T_true"), (1, "T_alt1")]),
    ("Old admix frac (d)", "d_frac", F_OLD,         [(0, "T_true"), (1, "T_alt1")]),
]

fig2, axes2 = plt.subplots(2, 3, figsize=(24, 12))

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
    f"Old-admix parameter recovery (v4)  |  "
    f"truth: t={T_OLD}, f={F_OLD}  |  Ne_d={NE_D}",
    fontsize=13, fontweight="bold"
)
fig2.tight_layout(rect=[0, 0, 1, 0.93])
fig2.savefig("mixed_beats_both_v4_params.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: mixed_beats_both_v4_params.png")


# ====================================================================
# 6. Summary table
# ====================================================================
print("\n" + "=" * 110)
print("MIXED BEATS BOTH v4 - SUMMARY")
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
