"""
Mixed-beats-both experiment v2:
  Weak admixture fractions (0.1) + varying Ne + extended IBD bins.

True topology (T_true): 5 populations (a, b, c, d, e) with TWO admixture events:
  - RECENT admixture on b (t=20): 10% from a-side (bP1), 90% from c-side (bP2)
  - OLD admixture on d (t=500): 10% from c-side (dP1), 90% from e-side (dP2)

  Going backward in time:
    t=20:   b splits -> bP1, bP2               (recent admixture)
    t=50:   a + bP1 -> ab
    t=100:  c + bP2 -> cb
    t=500:  d splits -> dP1, dP2               (old admixture)
    t=700:  cb + dP1 -> cbd
    t=1000: ab + cbd -> left
    t=1200: e + dP2 -> right
    t=1500: left + right -> root

  Ne: b=20000, d=5000, all others=10000 (haploid)

Alternative topologies:
  T_alt1 (old-only): only the OLD d-admixture, b is a normal sister of a
  T_alt2 (recent-only): only the RECENT b-admixture, d is a normal sister of c
"""

import sys, os
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
    """Extract the best ELBO from pathfinder stdout files."""
    best_elbo = -np.inf
    for stdout_file in fit._runset.stdout_files:
        try:
            with open(stdout_file, 'r') as f:
                for line in f:
                    m = re.search(r'Best Iter:.*?ELBO \(([-+\d.eE]+)\)', line)
                    if m:
                        elbo_val = float(m.group(1))
                        if elbo_val > best_elbo:
                            best_elbo = elbo_val
        except (OSError, ValueError):
            continue
    return best_elbo if np.isfinite(best_elbo) else np.nan


def get_admix_event_mapping(dem):
    """Return list of (child_name, cum_time_idx, admix_frac_idx) for each admixture event."""
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

# Extended IBD bins: 0.2 cM lower bound (tree-sequence MRCA can resolve this)
# Old admixture at t=500 produces segments ~0.1 cM on average; the tail above
# 0.2 cM is still informative.
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

# Population sizes (haploid)
NE_B = 5000
NE_D = 5000
NE_REST = 15000

# Admixture fractions: fraction_parent_1 = fraction to FIRST named parent
# b: 10% to bP1 (a-side), 90% to bP2 (c-side)
# d: 10% to dP1 (c-side), 90% to dP2 (e-side)
ADMIX_FRAC_RECENT = 0.9
ADMIX_FRAC_OLD = 0.7

SNP_CUTOFF_TIME = 1500
SNP_MIN_MAF = 0.05
FIXED_NE = NE_REST


# ====================================================================
# 1. Define topologies
# ====================================================================
print("=" * 60)
print("Defining topologies")
print("=" * 60)

# --- T_true: both admixture events ---
dem_sim = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
dem_sim.add_admixture_event('b', 'bP1', 'bP2')
dem_sim.add_merge_event('a', 'bP1', 'ab')
dem_sim.add_merge_event('c', 'bP2', 'cb')
dem_sim.add_admixture_event('d', 'dP1', 'dP2')
dem_sim.add_merge_event('cb', 'dP1', 'cbd')
dem_sim.add_merge_event('ab', 'cbd', 'left')
dem_sim.add_merge_event('e', 'dP2', 'right')
dem_sim.add_merge_event('left', 'right', 'root')

# Set Ne (diploid for msprime = haploid / 2)
for node in dem_sim.nodes:
    dem_sim.set_node_ne(node, NE_REST // 2)
dem_sim.set_node_ne('b', NE_B // 2)
dem_sim.set_node_ne('d', NE_D // 2)

# Admixture parameters
dem_sim.set_admixture_parameters('b', time=20,
                                 fraction_parent_1=ADMIX_FRAC_RECENT,
                                 parent_1_name='bP1')
dem_sim.set_admixture_parameters('d', time=500,
                                 fraction_parent_1=ADMIX_FRAC_OLD,
                                 parent_1_name='dP1')

# Merge times
dem_sim.set_merge_time('ab', 50)
dem_sim.set_merge_time('cb', 100)
dem_sim.set_merge_time('cbd', 700)
dem_sim.set_merge_time('left', 1000)
dem_sim.set_merge_time('right', 1200)
dem_sim.set_merge_time('root', 1500)
dem_sim.finalize_root()

dem_true = dem_sim

# --- T_alt1: old admixture only (no recent b-admixture) ---
# b is a normal sister of a; a+b → ab → ab+c → abc; then d admixes
dem_alt1 = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
dem_alt1.add_merge_event('a', 'b', 'ab')
dem_alt1.add_merge_event('ab', 'c', 'abc')
dem_alt1.add_admixture_event('d', 'dP1', 'dP2')
dem_alt1.add_merge_event('abc', 'dP1', 'left')
dem_alt1.add_merge_event('e', 'dP2', 'right')
dem_alt1.add_merge_event('left', 'right', 'root')
dem_alt1.finalize_root()

# --- T_alt2: recent admixture only (no old d-admixture) ---
dem_alt2 = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
dem_alt2.add_admixture_event('b', 'bP1', 'bP2')
dem_alt2.add_merge_event('a', 'bP1', 'ab')
dem_alt2.add_merge_event('c', 'bP2', 'cb')
dem_alt2.add_merge_event('cb', 'd', 'cbd')
dem_alt2.add_merge_event('ab', 'cbd', 'left')
dem_alt2.add_merge_event('left', 'e', 'root')
dem_alt2.finalize_root()

topology_dems = [dem_true, dem_alt1, dem_alt2]
topology_labels = [
    "T_true (both admix)",
    "T_alt1 (old only)",
    "T_alt2 (recent only)",
]

for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, n_admix={d.n_admix}, "
          f"n_nodes={len(d.nodes)}")

n_topos = len(topology_dems)


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
    mts = simulate_msprime(
        dem_sim,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )

    packed, n_blocks, pop_samples, pids, pi = calculate_ibd_blocks_mrca(
        mts,
        bins=bins,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
    )
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi
        pop_ids = pids

    snp_blks, _ = prepare_snp_blocks(
        mts, dem_sim,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME,
        min_maf=SNP_MIN_MAF,
    )
    snp_blocks_list.append(snp_blks)

pooled_ibd, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
pooled_snp = pool_snp_blocks(snp_blocks_list)

print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} independent simulations")


# ====================================================================
# 3. Compile Stan models (Nvarying)
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
for ti, (lbl, mapping) in enumerate(zip(topology_labels, admix_mappings)):
    print(f"  {lbl} admix mapping: {mapping}")

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
            pooled_ibd,
            n_blocks_total=total_blocks,
            pop_ids=pop_ids,
            pair_info=pair_info,
            bins=bins,
            target_cm=cm_val,
            block_size_cm=BLOCK_SIZE_CM,
            chosen_blocks=chosen_blocks,
        )

        leaves = dem_sim.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid,
            se_block_size=50,
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

                    # Extract admixture parameter posterior means
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
    ("T_true - T_alt1 (detects recent admix?)", 0, 1),
    ("T_true - T_alt2 (detects old admix?)", 0, 2),
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
    "Mixed beats both v2: T_true (2 admix) vs T_alt1 (old only) vs T_alt2 (recent only)\n"
    f"Ne: b={NE_B}, d={NE_D}, rest={NE_REST}  |  "
    f"admix fracs: recent={ADMIX_FRAC_RECENT}, old={ADMIX_FRAC_OLD}  |  "
    f"bins from {bins[0][0]} cM",
    fontsize=12, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("mixed_beats_both_v2.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: mixed_beats_both_v2.png")


# ====================================================================
# 5b. Plot: Parameter estimation (4x3 grid)
#     Rows: recent_time(b), recent_frac(b), old_time(d), old_frac(d)
#     Cols: IBD-only, SNP-only, Mixed
#     Colors: T_true = blue, T_alt1 = orange, T_alt2 = green
# ====================================================================
TOPO_COLORS = {
    0: "#1f77b4",   # T_true = blue
    1: "#ff7f0e",   # T_alt1 = orange
    2: "#2ca02c",   # T_alt2 = green
}
TOPO_NAMES = {0: "T_true", 1: "T_alt1", 2: "T_alt2"}

# (row_label, param_key, truth, [(topo_idx, ...)])
param_rows = [
    ("Recent admix time (b)", "b_time", 20.0, [
        (0, "T_true (both admix)"),
        (2, "T_alt2 (recent only)"),
    ]),
    ("Recent admix frac (b)", "b_frac", ADMIX_FRAC_RECENT, [
        (0, "T_true (both admix)"),
        (2, "T_alt2 (recent only)"),
    ]),
    ("Old admix time (d)", "d_time", 500.0, [
        (0, "T_true (both admix)"),
        (1, "T_alt1 (old only)"),
    ]),
    ("Old admix frac (d)", "d_frac", ADMIX_FRAC_OLD, [
        (0, "T_true (both admix)"),
        (1, "T_alt1 (old only)"),
    ]),
]

fig2, axes2 = plt.subplots(4, 3, figsize=(24, 22))

for row_idx, (row_label, param_key, truth, topo_entries) in enumerate(param_rows):
    for m_idx, m_label in enumerate(model_labels):
        ax = axes2[row_idx, m_idx]
        n_g = len(topo_entries)
        bar_width = 0.7 / n_g

        for gi, (ti, topo_disp) in enumerate(topo_entries):
            color = TOPO_COLORS[ti]
            data_list = []
            positions = []
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

        # Build descriptive title showing which topologies are compared
        topo_desc = " vs ".join(
            f"{TOPO_NAMES[ti]}" for ti, _ in topo_entries
        )
        ax.set_title(f"{m_label}: {row_label}\n({topo_desc})",
                     fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.15)

        # Legend on first column only
        if m_idx == 0:
            handles = [Patch(facecolor=TOPO_COLORS[ti], alpha=0.55,
                             label=topo_disp) for ti, topo_disp in topo_entries]
            handles.append(plt.Line2D([0], [0], color='red', ls='--',
                                       lw=1.2, label="truth"))
            ax.legend(handles=handles, fontsize=7, loc='best')

fig2.suptitle(
    "Admixture parameter estimation (v2)\n"
    f"True: b admix at t=20 (frac={ADMIX_FRAC_RECENT}), "
    f"d admix at t=500 (frac={ADMIX_FRAC_OLD})  |  "
    f"Ne: b={NE_B}, d={NE_D}, rest={NE_REST}",
    fontsize=13, fontweight="bold"
)
fig2.tight_layout(rect=[0, 0, 1, 0.94])
fig2.savefig("mixed_beats_both_v2_params.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: mixed_beats_both_v2_params.png")


# ====================================================================
# 6. Summary table
# ====================================================================
print("\n" + "=" * 110)
print("MIXED BEATS BOTH v2 — SUMMARY")
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
