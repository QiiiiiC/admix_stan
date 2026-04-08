"""
Run N_REPLICATES of IBD-based demographic inference using:
  - Multiple independent simulations pooled for block-bootstrap
  - Delete-one-haplotype jackknife for variance estimation
  - CmdStan MCMC sampling (full posterior, no pathfinder approximation)

Purpose: diagnose whether the ~2.5% downward Ne bias seen with pathfinder
is an inference artifact or a model issue.
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from simulation_methods import simulate_msprime, build_ibd_stan_data
from ibd_jackknife import (
    calculate_ibd_blocks_mrca,
    pool_multiple_simulations,
    resample_ibd_with_jackknife_variance,
)
from demography import DemographicTopology
from cmdstanpy import CmdStanModel


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 50  # fewer replicates since MCMC is slower
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}

# Number of independent simulations to pool
N_SIMS = 5
SIM_CM_EACH = 500  # each sim is 500 cM -> 5 x 500 = 2500 cM total, 50 blocks
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 200, 500]  # fewer cm values for speed

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 10.0], [10.0, max(cm_values)]
]

true_vals = {
    "Admixture time": 20,
    "Effective population size": 6000,
    "Admixture fraction": 0.75,
    "Post-admixture merge time": 200,
}

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
dem.set_merge_time('subAnc', 200)
dem.set_merge_time('anc', 500)
dem.finalize_root()

# ====================================================================
# 2. Simulate MULTIPLE independent genomes and pool blocks
# ====================================================================
packed_list = []
n_blocks_list = []
pair_info = None

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    mts = simulate_msprime(
        dem,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,  # different seed for each
    )

    packed, n_blocks, pop_samples, pop_ids, pi = calculate_ibd_blocks_mrca(
        mts,
        bins=bins,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
    )
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi  # same for all sims (same topology & sample sizes)

# Pool all blocks into one array
pooled, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} independent simulations")

# ====================================================================
# 3. Compile Stan model once
# ====================================================================
model = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nfixed.stan"))

# ====================================================================
# 4. Run replicates with MCMC
# ====================================================================
results_mcmc = {cm_val: [] for cm_val in cm_values}
results_pathfinder = {cm_val: [] for cm_val in cm_values}
rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} replicates)")
    print(f"{'='*60}")

    for rep in range(N_REPLICATES):
        if (rep + 1) % 5 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        ibd_mean, ibd_var = resample_ibd_with_jackknife_variance(
            pooled,
            n_blocks_total=total_blocks,
            pop_ids=pop_ids,
            pair_info=pair_info,
            bins=bins,
            target_cm=cm_val,
            block_size_cm=BLOCK_SIZE_CM,
            rng=rng,
        )

        stan_ibd = build_ibd_stan_data(
            dem, ibd_mean, ibd_var, bins,
            T_max=100000, cm=cm_val,
        )

        init_dict = {
            "times": [100.0] * stan_ibd['n_events'],
            "admixture_fractions": [0.5] * stan_ibd['n_admixture'],
            "effective_N": 10000.0,
        }

        # --- MCMC sampling ---
        try:
            fit_mcmc = model.sample(
                data=stan_ibd,
                inits=init_dict,
                chains=2,
                iter_warmup=500,
                iter_sampling=500,
                show_console=False,
                seed=rep + cm_val,
            )

            all_vars = fit_mcmc.stan_variables()
            pmeans = {
                name: draws.mean(axis=0)
                for name, draws in all_vars.items()
            }
            results_mcmc[cm_val].append(pmeans)

        except Exception as e:
            print(f"    [WARN] MCMC replicate {rep+1} failed: {e}")

        # --- Pathfinder (same data, for comparison) ---
        try:
            fit_pf = model.pathfinder(
                data=stan_ibd,
                inits=init_dict,
                show_console=False,
            )

            all_vars = fit_pf.stan_variables()
            pmeans = {
                name: draws.mean(axis=0)
                for name, draws in all_vars.items()
            }
            results_pathfinder[cm_val].append(pmeans)

        except Exception as e:
            print(f"    [WARN] Pathfinder replicate {rep+1} failed: {e}")

# ====================================================================
# 5. Extract posterior means for plotting
# ====================================================================
def extract_param(pmeans_dict, param_name):
    if param_name == "Admixture time":
        return pmeans_dict['times'][0] if np.ndim(pmeans_dict['times']) > 0 else pmeans_dict['times']
    elif param_name == "Effective population size":
        return float(pmeans_dict['effective_N'])
    elif param_name == "Admixture fraction":
        af = pmeans_dict['admixture_fractions']
        return af[0] if np.ndim(af) > 0 else float(af)
    elif param_name == "Post-admixture merge time":
        t = pmeans_dict['times']
        return float(np.sum(t[:4]))
    else:
        raise ValueError(f"Unknown param: {param_name}")

# ====================================================================
# 6. Plot side-by-side comparison
# ====================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, (name, tv) in zip(axes.flatten(), true_vals.items()):
    mcmc_samples = []
    pf_samples = []

    for cm_val in cm_values:
        mcmc_vals = [extract_param(pm, name) for pm in results_mcmc[cm_val]]
        pf_vals = [extract_param(pm, name) for pm in results_pathfinder[cm_val]]
        mcmc_samples.append(np.array(mcmc_vals))
        pf_samples.append(np.array(pf_vals))

    positions = np.arange(len(cm_values))
    width = 0.3

    # MCMC boxplots (left, blue)
    bp1 = ax.boxplot(
        mcmc_samples,
        positions=positions - width/2,
        widths=width * 0.8,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="#BA7517", linewidth=2),
        whiskerprops=dict(color="k", linewidth=1.2),
        capprops=dict(color="k", linewidth=1.2),
        boxprops=dict(facecolor="#378ADD", alpha=0.4, edgecolor="#185FA5", linewidth=1.2),
    )

    # Pathfinder boxplots (right, orange)
    bp2 = ax.boxplot(
        pf_samples,
        positions=positions + width/2,
        widths=width * 0.8,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="#BA7517", linewidth=2),
        whiskerprops=dict(color="k", linewidth=1.2),
        capprops=dict(color="k", linewidth=1.2),
        boxprops=dict(facecolor="#E8A838", alpha=0.4, edgecolor="#C07010", linewidth=1.2),
    )

    ax.axhline(tv, color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}")
    ax.legend(
        [bp1["boxes"][0], bp2["boxes"][0], ax.get_lines()[-1]],
        ["MCMC", "Pathfinder", f"True = {tv}"],
        fontsize=9,
    )
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{cm} cM" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_title(name)
    ax.grid(axis="y", alpha=0.15)

plt.suptitle("MCMC vs Pathfinder: Ne bias diagnostic", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("mcmc_vs_pathfinder_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: mcmc_vs_pathfinder_comparison.png")

# Print summary statistics
print("\n" + "="*60)
print("Summary: Mean Ne across replicates")
print("="*60)
for cm_val in cm_values:
    mcmc_ne = [extract_param(pm, "Effective population size") for pm in results_mcmc[cm_val]]
    pf_ne = [extract_param(pm, "Effective population size") for pm in results_pathfinder[cm_val]]
    print(f"  {cm_val:4d} cM: MCMC = {np.mean(mcmc_ne):8.1f} (n={len(mcmc_ne)}), "
          f"Pathfinder = {np.mean(pf_ne):8.1f} (n={len(pf_ne)}), "
          f"True = 6000")
