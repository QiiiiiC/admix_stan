"""
Compare three models using synchronized block bootstrap:
  1. IBD-only   (ibd_model_Nfixed.stan, effective_N as parameter)
  2. SNP-only   (snp_model_Nfixed.stan, effective_N fixed at 6000)
  3. Mixed       (mixed_model.stan, effective_N as parameter)

Same genomic blocks are selected for both IBD and SNP resampling.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from simulation_methods import simulate_msprime, build_ibd_stan_data, build_snp_stan_data, build_mixed_stan_data
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
import re


def extract_elbo(fit):
    """Extract the best ELBO from pathfinder stdout files.

    CmdStan prints: Path [i] : Best Iter: [N] ELBO (value) evaluations: (count)
    This parses those lines and returns the maximum ELBO across paths.
    Works regardless of psis_resample setting.
    """
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

# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 200
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}

# Number of independent simulations to pool
N_SIMS = 50
SIM_CM_EACH = 50  # each sim is 50 cM → 50 x 50 = 2500 cM total, 50 blocks
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65],[0.65, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],[8.0,20.0], [20.0, BLOCK_SIZE_CM]
]

true_vals = {
    "Admixture time": 20,
    "Effective population size": 6000,
    "Admixture fraction": 0.75,
    "Post-admixture merge time": 1000,
}

SNP_CUTOFF_TIME = 1500  # keep only mutations older than this
SNP_MIN_MAF = 0.05
FIXED_NE = 6000         # fixed Ne for SNP-only model

# ====================================================================
# 1. Define topology
# ====================================================================
dem = DemographicTopology(['a', 'admix', 'b', 'c'])
dem.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem.add_merge_event('a', 'admixPa1', 'aAdmix')
dem.add_merge_event('b', 'admixPa2', 'bAdmix')
dem.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
dem.add_merge_event('subAnc', 'c', 'anc')
for node in ['a','b','c','admix','admixPa1','admixPa2','aAdmix','bAdmix','subAnc','anc']:
    dem.set_node_ne(node, 3000)
dem.set_admixture_parameters('admix', 20, 0.75, 'admixPa1')
dem.set_merge_time('aAdmix', 40)
dem.set_merge_time('bAdmix', 50)
dem.set_merge_time('subAnc', 1000)
dem.set_merge_time('anc', 1500)
dem.finalize_root()

# ====================================================================
# 2. Simulate and pool blocks (same sims for both IBD and SNP)
# ====================================================================
packed_list = []
n_blocks_list = []
snp_blocks_list = []
pair_info = None

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    mts = simulate_msprime(
        dem,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )

    # IBD blocks
    packed, n_blocks, pop_samples, pop_ids, pi = calculate_ibd_blocks_mrca(
        mts,
        bins=bins,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
    )
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi

    # SNP blocks: per-SNP (deviation, heterozygosity) grouped by genomic block
    snp_blks, _ = prepare_snp_blocks(
        mts, dem,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME,
        min_maf=SNP_MIN_MAF,
    )
    snp_blocks_list.append(snp_blks)

# Pool across simulations
pooled_ibd, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
pooled_snp = pool_snp_blocks(snp_blocks_list)

print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} independent simulations")

# ====================================================================
# 3. Compile Stan models
# ====================================================================
ibd_model = CmdStanModel(stan_file="ibd_model_Nfixed.stan")
snp_model = CmdStanModel(stan_file="snp_model_Nfixed.stan")
mixed_model = CmdStanModel(stan_file="mixed_model.stan")

# ====================================================================
# 4. Run replicates
# ====================================================================
results_ibd = {cm_val: [] for cm_val in cm_values}
results_snp = {cm_val: [] for cm_val in cm_values}
results_mixed = {cm_val: [] for cm_val in cm_values}
elbo_ibd = {cm_val: [] for cm_val in cm_values}
elbo_snp = {cm_val: [] for cm_val in cm_values}
elbo_mixed = {cm_val: [] for cm_val in cm_values}

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} replicates)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        # Select blocks ONCE for both IBD and SNP (without replacement)
        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed, replace=False)

        # --- Resample IBD ---
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

        # --- Resample SNP (TreeMix: gather SNPs, 50-SNP SE blocks) ---
        leaves = dem.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid,
            se_block_size=50,
        )

        # --- Init dicts ---
        n_events = len(dem.ordered_events)
        n_admix = dem.n_admix

        init_ibd = {
            "times": [100.0] * n_events,
            "admixture_fractions": [0.5] * n_admix,
            "effective_N": 10000.0,
        }
        init_snp = {
            "times": [100.0] * n_events,
            "admixture_fractions": [0.5] * n_admix,
        }
        init_mixed = {
            "times": [100.0] * n_events,
            "admixture_fractions": [0.5] * n_admix,
            "effective_N": 10000.0,
            "kappa_snp": 1.0,
        }

        # --- Model 1: IBD-only ---
        try:
            stan_data = build_ibd_stan_data(
                dem, ibd_mean, ibd_var, bins,
                T_max=100000, cm=cm_val,
            )
            fit = ibd_model.pathfinder(
                data=stan_data, inits=init_ibd, show_console=False,
                psis_resample=True,
            )
            all_vars = fit.stan_variables()
            pmeans = {name: draws.mean(axis=0) for name, draws in all_vars.items()}
            results_ibd[cm_val].append(pmeans)
            elbo_ibd[cm_val].append(extract_elbo(fit))
        except Exception as e:
            print(f"    [WARN] IBD rep {rep+1} failed: {e}")

        # --- Model 2: SNP-only ---
        try:
            stan_data = build_snp_stan_data(
                dem, w_hat, w_se, effective_N=FIXED_NE,
            )
            fit = snp_model.pathfinder(
                data=stan_data, inits=init_snp, show_console=False,
                psis_resample=True,
            )
            all_vars = fit.stan_variables()
            pmeans = {name: draws.mean(axis=0) for name, draws in all_vars.items()}
            results_snp[cm_val].append(pmeans)
            elbo_snp[cm_val].append(extract_elbo(fit))
        except Exception as e:
            print(f"    [WARN] SNP rep {rep+1} failed: {e}")

        # --- Model 3: Mixed ---
        try:
            stan_data = build_mixed_stan_data(
                dem, ibd_mean, ibd_var, bins, w_hat, w_se,
                T_max=100000, cm=cm_val,
            )
            fit = mixed_model.pathfinder(
                data=stan_data, inits=init_mixed, show_console=False,
                psis_resample=True,
            )
            all_vars = fit.stan_variables()
            pmeans = {name: draws.mean(axis=0) for name, draws in all_vars.items()}
            results_mixed[cm_val].append(pmeans)
            elbo_mixed[cm_val].append(extract_elbo(fit))
        except Exception as e:
            print(f"    [WARN] Mixed rep {rep+1} failed: {e}")


# ====================================================================
# 5. Extract parameters
# ====================================================================
def extract_param(pmeans_dict, param_name, has_ne=True):
    if param_name == "Admixture time":
        t = pmeans_dict['times']
        return t[0] if np.ndim(t) > 0 else float(t)
    elif param_name == "Effective population size":
        if has_ne:
            return float(pmeans_dict['effective_N'])
        else:
            return np.nan
    elif param_name == "Admixture fraction":
        af = pmeans_dict['admixture_fractions']
        return af[0] if np.ndim(af) > 0 else float(af)
    elif param_name == "Post-admixture merge time":
        t = pmeans_dict['times']
        return float(np.sum(t[:4]))
    else:
        raise ValueError(f"Unknown param: {param_name}")


# ====================================================================
# 6a. Plot 1: Parameter comparison (2x2)
# ====================================================================
fig1, axes1 = plt.subplots(2, 2, figsize=(16, 12))

model_labels = ["IBD-only", "SNP-only", "Mixed"]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]
model_results = [results_ibd, results_snp, results_mixed]
model_has_ne = [True, False, True]
from matplotlib.patches import Patch

for ax, (name, tv) in zip(axes1.flatten(), true_vals.items()):
    if name == "Effective population size":
        active_labels = ["IBD-only", "Mixed"]
        active_colors = ["#378ADD", "#E8A838"]
        active_results = [results_ibd, results_mixed]
        active_has_ne = [True, True]
    else:
        active_labels = model_labels
        active_colors = model_colors
        active_results = model_results
        active_has_ne = model_has_ne

    n_models = len(active_labels)
    n_cms = len(cm_values)

    positions_all = []
    data_all = []
    colors_all = []

    for cm_idx, cm_val in enumerate(cm_values):
        group_center = cm_idx * (n_models + 1)
        for m_idx in range(n_models):
            pos = group_center + m_idx
            positions_all.append(pos)
            colors_all.append(active_colors[m_idx])

            vals = [extract_param(pm, name, active_has_ne[m_idx])
                    for pm in active_results[m_idx][cm_val]]
            data_all.append(np.array(vals))

    bp = ax.boxplot(
        data_all,
        positions=positions_all,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        whiskerprops=dict(color="k", linewidth=1),
        capprops=dict(color="k", linewidth=1),
    )

    for patch, color in zip(bp['boxes'], colors_all):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    ax.axhline(tv, color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}")

    group_centers = [cm_idx * (n_models + 1) + (n_models - 1) / 2
                     for cm_idx in range(n_cms)]
    ax.set_xticks(group_centers)
    ax.set_xticklabels([f"{cm}" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_title(name)
    ax.grid(axis="y", alpha=0.15)

    handles = [Patch(facecolor=c, alpha=0.5, label=l)
               for c, l in zip(active_colors, active_labels)]
    handles.append(plt.Line2D([0], [0], color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}"))
    ax.legend(handles=handles, fontsize=8)

fig1.suptitle("IBD-only vs SNP-only vs Mixed: Posterior Mean Comparison",
              fontsize=14, fontweight="bold")
fig1.tight_layout()
fig1.savefig("snp_ibd_mixed_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: snp_ibd_mixed_comparison.png")

# ====================================================================
# 6b. Plot 2: ELBO comparison (overlapping boxplots per genome length)
# ====================================================================
fig2, ax_elbo = plt.subplots(1, 1, figsize=(14, 5))

elbo_results = [elbo_ibd, elbo_snp, elbo_mixed]
n_models = len(model_labels)
n_cms = len(cm_values)

# Plot each model as a separate boxplot layer at the same x positions
x_positions = list(range(n_cms))
for m_idx in range(n_models):
    data = [np.array(elbo_results[m_idx][cm_val]) for cm_val in cm_values]
    bp = ax_elbo.boxplot(
        data,
        positions=x_positions,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="k", linewidth=1.5),
        whiskerprops=dict(color=model_colors[m_idx], linewidth=1, alpha=0.6),
        capprops=dict(color=model_colors[m_idx], linewidth=1, alpha=0.6),
    )
    for patch in bp['boxes']:
        patch.set_facecolor(model_colors[m_idx])
        patch.set_alpha(0.3)

    # Trend line connecting mean ELBO
    mean_values = [np.mean(d) for d in data]
    ax_elbo.plot(x_positions, mean_values, color=model_colors[m_idx],
                 linewidth=2, marker='o', markersize=5, zorder=5,
                 label=model_labels[m_idx])

ax_elbo.set_xticks(x_positions)
ax_elbo.set_xticklabels([f"{cm}" for cm in cm_values])
ax_elbo.set_xlabel("Genome length (cM)")
ax_elbo.set_ylabel("ELBO")
ax_elbo.set_title("ELBO (from pathfinder optimization)")
ax_elbo.grid(axis="y", alpha=0.15)
ax_elbo.legend(fontsize=8)

fig2.tight_layout()
fig2.savefig("elbo_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: elbo_comparison.png")

# ====================================================================
# 7. Print summary table
# ====================================================================
print("\n" + "=" * 100)
print("SUMMARY TABLE — Bias (%) of posterior mean")
print("=" * 100)
print(f"{'cm':>6s} | {'Model':>10s} | {'Ne bias%':>10s} | {'Tadmix bias%':>12s} | "
      f"{'frac bias%':>10s} | {'Tmerge bias%':>12s} | {'n_ok':>5s}")
print("-" * 100)

for cm_val in cm_values:
    for label, res, has_ne in [
        ("IBD", results_ibd, True),
        ("SNP", results_snp, False),
        ("Mixed", results_mixed, True),
    ]:
        biases = {}
        for pname, tv in true_vals.items():
            vals = np.array([extract_param(pm, pname, has_ne) for pm in res[cm_val]])
            if pname == "Effective population size" and not has_ne:
                biases[pname] = np.nan
            else:
                biases[pname] = (np.nanmean(vals) - tv) / tv * 100

        ne_str = f"{biases['Effective population size']:+9.2f}%" if has_ne else "    fixed"
        print(f"{cm_val:6d} | {label:>10s} | {ne_str} | "
              f"{biases['Admixture time']:+11.2f}% | "
              f"{biases['Admixture fraction']:+9.2f}% | "
              f"{biases['Post-admixture merge time']:+11.2f}% | "
              f"{len(res[cm_val]):5d}")
