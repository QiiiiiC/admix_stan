"""
Model selection experiment: Can ELBO distinguish the correct demographic topology?

Data is simulated under T1 (true: admix = 75% a + 25% b). Three candidate
topologies are fit to the same resampled data using IBD-only, SNP-only,
and Mixed models.

Topologies:
  T1 (true):      admix is admixture of a(75%) and b(25%)
  T2 (wrong src): admix is admixture of a(75%) and c(25%)
  T3 (no admix):  admix is tree-like sister to a (no admixture)

Plots:
  1. ELBO(T1) - ELBO(T_wrong) vs genome length (3 panels: IBD, SNP, Mixed)
  2. Softmax model weight w(T1) vs genome length (1 panel, 3 model types)
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
import re
from contextlib import redirect_stdout

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
from relative_error import plot_relative_error_boxplot


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


def softmax_weights(elbo_array):
    """Compute softmax model weights from ELBO values (numerically stable)."""
    max_elbo = np.max(elbo_array, axis=1, keepdims=True)
    exp_elbo = np.exp(elbo_array - max_elbo)
    return exp_elbo / exp_elbo.sum(axis=1, keepdims=True)


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 200
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}
# Haploid sample count per population (diploid samples x ploidy 2)
N_HAPLOID_PER_POP = {_p: 2 * _c for _p, _c in SAMPLES_PER_POP.items()}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65],[0.65, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],[8.0,20.0], [20.0, BLOCK_SIZE_CM]
]

SNP_CUTOFF_TIME = 1500
SNP_MIN_MAF = 0.05


# ====================================================================
# 1. Define topologies
# ====================================================================
print("=" * 60)
print("Defining topologies")
print("=" * 60)

# --- Simulation topology (T1 with full parameters for msprime) ---
dem_sim = DemographicTopology(['a', 'admix', 'b', 'c'])
dem_sim.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem_sim.add_merge_event('a', 'admixPa1', 'aAdmix')
dem_sim.add_merge_event('b', 'admixPa2', 'bAdmix')
dem_sim.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
dem_sim.add_merge_event('subAnc', 'c', 'anc')
for node in dem_sim.nodes:
    dem_sim.set_node_ne(node, 3000)
dem_sim.set_admixture_parameters('admix', 20, 0.75, 'admixPa1')
dem_sim.set_merge_time('aAdmix', 40)
dem_sim.set_merge_time('bAdmix', 50)
dem_sim.set_merge_time('subAnc', 1000)
dem_sim.set_merge_time('anc', 1500)
dem_sim.finalize_root()

# T1 (true): same structure — reuse for inference
dem_T1 = dem_sim

# T2 (wrong source): admix from a and c instead of a and b
print("\nT2 (wrong source): admix from a and c")
dem_T2 = DemographicTopology(['a', 'admix', 'b', 'c'])
dem_T2.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem_T2.add_merge_event('a', 'admixPa1', 'aAdmix')
dem_T2.add_merge_event('c', 'admixPa2', 'cAdmix')
dem_T2.add_merge_event('aAdmix', 'b', 'subAnc')
dem_T2.add_merge_event('subAnc', 'cAdmix', 'anc')
dem_T2.finalize_root()

# T3 (no admixture): admix is tree-like sister to a
print("\nT3 (no admixture): admix is sister to a")
dem_T3 = DemographicTopology(['a', 'admix', 'b', 'c'])
dem_T3.add_merge_event('admix', 'a', 'admixA')
dem_T3.add_merge_event('admixA', 'b', 'admixAB')
dem_T3.add_merge_event('admixAB', 'c', 'anc')
dem_T3.finalize_root()

topology_dems = [dem_T1, dem_T2, dem_T3]
topology_labels = ["T1 (true)", "T2 (wrong src)", "T3 (no admix)"]
n_topologies = len(topology_dems)

for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, n_admix={d.n_admix}, "
          f"n_nodes={len(d.nodes)}")


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
# 3. Compile Stan models
# ====================================================================
ibd_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nfixed.stan"))
snp_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "snp_model_Nfixed.stan"))
mixed_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nfixed.stan"))


# ====================================================================
# 4. Run replicates: fit 3 topologies × 3 model types
# ====================================================================
model_types = ["ibd", "snp", "mixed"]
model_labels_type = ["IBD-only", "SNP-only", "Mixed"]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]

# Aligned ELBO storage: elbo_all[model_type][cm_val] = list of (e_T1, e_T2, e_T3)
# Only stores replicates where ALL 3 topologies succeeded for that model type
elbo_all = {mt: {cm: [] for cm in cm_values} for mt in model_types}

# Relative-error collector: Mixed model, true topology (ti==0), largest cm only
rel_err_stanvars = []
REL_ERR_CM = cm_values[-1]

rng = np.random.default_rng(seed=2025)
_devnull = open(os.devnull, 'w')  # suppress verbose build output

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} reps × {n_topologies} topologies × 3 models)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        # Select blocks ONCE (shared across all models and topologies)
        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed, replace=False)

        # --- Resample data (shared) ---
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

        # --- Fit each topology, collect ELBO per model type ---
        rep_elbos = {mt: [np.nan] * n_topologies for mt in model_types}

        for ti, dem_infer in enumerate(topology_dems):
            n_events = len(dem_infer.ordered_events)
            n_admix = dem_infer.n_admix

            init_base = {"times": [100.0] * n_events}
            if n_admix > 0:
                init_base["admixture_fractions"] = [0.5] * n_admix

            # --- IBD-only ---
            try:
                with redirect_stdout(_devnull):
                    sd = build_ibd_stan_data(
                        dem_infer, ibd_mean, ibd_var, bins, n_samples_per_pop=N_HAPLOID_PER_POP,
                        T_max=100000, cm=cm_val,
                    )
                init = {**init_base, "effective_N": 10000.0}
                fit = ibd_stan.pathfinder(
                    data=sd, inits=init, show_console=False, psis_resample=True,
                )
                rep_elbos["ibd"][ti] = extract_elbo(fit)
            except Exception as e:
                if rep < 3:
                    print(f"    [WARN] IBD {topology_labels[ti]} rep {rep+1}: {e}")

            # --- SNP-only ---
            try:
                with redirect_stdout(_devnull):
                    sd = build_snp_stan_data(
                        dem_infer, w_hat, w_se,
                    )
                # SNP-only Nfixed has effective_N as a parameter; init it near the
                # truth or Pathfinder fails from its random start.
                init = {**init_base, "effective_N": 10000.0}
                fit = snp_stan.pathfinder(
                    data=sd, inits=init, show_console=False, psis_resample=True,
                )
                rep_elbos["snp"][ti] = extract_elbo(fit)
            except Exception as e:
                if rep < 3:
                    print(f"    [WARN] SNP {topology_labels[ti]} rep {rep+1}: {e}")

            # --- Mixed ---
            try:
                with redirect_stdout(_devnull):
                    sd = build_mixed_stan_data(
                        dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se, n_samples_per_pop=N_HAPLOID_PER_POP,
                        T_max=100000, cm=cm_val,
                    )
                init = {**init_base, "effective_N": 10000.0, "kappa_snp": 1.0}
                fit = mixed_stan.pathfinder(
                    data=sd, inits=init, show_console=False, psis_resample=True,
                )
                rep_elbos["mixed"][ti] = extract_elbo(fit)
                # Collect for relative-error plot: Mixed + true topology + largest cm
                if ti == 0 and cm_val == REL_ERR_CM:
                    rel_err_stanvars.append(fit.stan_variables())
            except Exception as e:
                if rep < 3:
                    print(f"    [WARN] Mixed {topology_labels[ti]} rep {rep+1}: {e}")

        # Store only replicates where ALL topologies succeeded for a model type
        for mt in model_types:
            if all(np.isfinite(rep_elbos[mt])):
                elbo_all[mt][cm_val].append(tuple(rep_elbos[mt]))

_devnull.close()


# ====================================================================
# 5. Plot 1: ELBO differences (T_true - T_wrong)
# ====================================================================
fig1, axes1 = plt.subplots(1, 3, figsize=(18, 5), sharey=False)

diff_colors = ["#378ADD", "#E24B4A"]

for ax, mt, mt_label in zip(axes1, model_types, model_labels_type):
    diff_T2_data = []
    diff_T3_data = []

    for cm_val in cm_values:
        elbos = elbo_all[mt][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)
            diff_T2_data.append(arr[:, 0] - arr[:, 1])
            diff_T3_data.append(arr[:, 0] - arr[:, 2])
        else:
            diff_T2_data.append(np.array([np.nan]))
            diff_T3_data.append(np.array([np.nan]))

    x = np.arange(len(cm_values))

    bp1 = ax.boxplot(
        diff_T2_data,
        positions=x - 0.18,
        widths=0.3,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        whiskerprops=dict(color=diff_colors[0], linewidth=1),
        capprops=dict(color=diff_colors[0], linewidth=1),
    )
    for patch in bp1['boxes']:
        patch.set_facecolor(diff_colors[0])
        patch.set_alpha(0.4)

    bp2 = ax.boxplot(
        diff_T3_data,
        positions=x + 0.18,
        widths=0.3,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        whiskerprops=dict(color=diff_colors[1], linewidth=1),
        capprops=dict(color=diff_colors[1], linewidth=1),
    )
    for patch in bp2['boxes']:
        patch.set_facecolor(diff_colors[1])
        patch.set_alpha(0.4)

    # Trend lines connecting means
    mean_T2 = [np.nanmean(d) for d in diff_T2_data]
    mean_T3 = [np.nanmean(d) for d in diff_T3_data]
    ax.plot(x, mean_T2, color=diff_colors[0], linewidth=2, marker='o', markersize=4,
            zorder=5, label="T1 \u2212 T2 (wrong src)")
    ax.plot(x, mean_T3, color=diff_colors[1], linewidth=2, marker='s', markersize=4,
            zorder=5, label="T1 \u2212 T3 (no admix)")

    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{cm}" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_ylabel("\u0394ELBO")
    ax.set_title(mt_label)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.15)

fig1.suptitle("ELBO(T_true) \u2212 ELBO(T_wrong): positive = correct topology preferred",
              fontsize=12, fontweight="bold")
fig1.tight_layout()
fig1.savefig("model_selection_elbo_difference.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: model_selection_elbo_difference.png")


# ====================================================================
# 6. Plot 2: Model weights w(T_true) — softmax of ELBO
# ====================================================================
fig2, ax2 = plt.subplots(1, 1, figsize=(14, 5))

x = np.arange(len(cm_values))

for mt_idx, (mt, mt_label, color) in enumerate(
    zip(model_types, model_labels_type, model_colors)
):
    weights_data = []
    for cm_val in cm_values:
        elbos = elbo_all[mt][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)             # (n_reps, 3)
            w = softmax_weights(arr)           # (n_reps, 3)
            weights_data.append(w[:, 0])       # w(T1)
        else:
            weights_data.append(np.array([1.0 / n_topologies]))

    bp = ax2.boxplot(
        weights_data,
        positions=x,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        whiskerprops=dict(color=color, linewidth=1, alpha=0.6),
        capprops=dict(color=color, linewidth=1, alpha=0.6),
    )
    for patch in bp['boxes']:
        patch.set_facecolor(color)
        patch.set_alpha(0.3)

    # Trend line
    mean_w = [np.mean(d) for d in weights_data]
    ax2.plot(x, mean_w, color=color, linewidth=2, marker='o', markersize=5,
             zorder=5, label=mt_label)

ax2.axhline(1.0 / n_topologies, color='k', ls='--', lw=1, alpha=0.5,
            label=f"Chance ({1.0/n_topologies:.2f})")
ax2.set_xticks(x)
ax2.set_xticklabels([f"{cm}" for cm in cm_values])
ax2.set_xlabel("Genome length (cM)")
ax2.set_ylabel("w(T_true)")
ax2.set_title("Model weight for true topology (softmax of ELBO)")
ax2.set_ylim(-0.05, 1.05)
ax2.legend(fontsize=9)
ax2.grid(axis="y", alpha=0.15)

fig2.tight_layout()
fig2.savefig("model_selection_weights.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: model_selection_weights.png")

# Relative-error box plot: Mixed model, true topology (T1), largest genome length
plot_relative_error_boxplot(
    dem_T1, rel_err_stanvars, varying=False,
    outpath="model_selection_relative_error_Nfixed.png",
    title="Parameter relative error (true topology)",
)


# ====================================================================
# 7. Summary table
# ====================================================================
print("\n" + "=" * 90)
print("MODEL SELECTION SUMMARY")
print("=" * 90)
print(f"{'cm':>6s} | {'Model':>10s} | {'w(T1)':>8s} | {'dE(T1-T2)':>12s} | "
      f"{'dE(T1-T3)':>12s} | {'n_ok':>5s}")
print("-" * 90)

for cm_val in cm_values:
    for mt, mt_label in zip(model_types, model_labels_type):
        elbos = elbo_all[mt][cm_val]
        if len(elbos) > 0:
            arr = np.array(elbos)
            w = softmax_weights(arr)
            mean_w = np.mean(w[:, 0])
            mean_diff_T2 = np.mean(arr[:, 0] - arr[:, 1])
            mean_diff_T3 = np.mean(arr[:, 0] - arr[:, 2])
            print(f"{cm_val:6d} | {mt_label:>10s} | {mean_w:8.3f} | "
                  f"{mean_diff_T2:+11.1f} | {mean_diff_T3:+11.1f} | {len(elbos):5d}")
        else:
            print(f"{cm_val:6d} | {mt_label:>10s} | {'N/A':>8s} | "
                  f"{'N/A':>12s} | {'N/A':>12s} | {0:5d}")
