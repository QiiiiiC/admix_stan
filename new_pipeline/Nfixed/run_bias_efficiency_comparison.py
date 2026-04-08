"""
Compare variational inference (pathfinder) vs MCMC:
  - Bias of posterior mean (% of true value) for 4 parameters
  - Running time

100 replicates per genome length, 5 genome lengths.
MCMC: chains=4, iter_warmup=1000, iter_sampling=2000.
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import time
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
N_REPLICATES = 100
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'admix': 15}

N_SIMS = 5
SIM_CM_EACH = 500
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 200, 300, 500]

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
# 2. Simulate and pool blocks
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
        seed=42 + sim_i,
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
        pair_info = pi

pooled, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM "
      f"= {total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} independent simulations")

# ====================================================================
# 3. Compile Stan model
# ====================================================================
model = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nfixed.stan"))

# ====================================================================
# 4. Run replicates
# ====================================================================
results_mcmc = {cm_val: [] for cm_val in cm_values}
results_pf = {cm_val: [] for cm_val in cm_values}
times_mcmc = {cm_val: [] for cm_val in cm_values}
times_pf = {cm_val: [] for cm_val in cm_values}

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} replicates)")
    print(f"{'='*60}")

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
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

        # --- MCMC ---
        try:
            t0 = time.time()
            fit_mcmc = model.sample(
                data=stan_ibd,
                inits=init_dict,
                chains=4,
                iter_warmup=1000,
                iter_sampling=2000,
                show_console=False,
                seed=rep + cm_val,
            )
            elapsed_mcmc = time.time() - t0

            all_vars = fit_mcmc.stan_variables()
            pmeans = {
                name: draws.mean(axis=0)
                for name, draws in all_vars.items()
            }
            results_mcmc[cm_val].append(pmeans)
            times_mcmc[cm_val].append(elapsed_mcmc)

        except Exception as e:
            print(f"    [WARN] MCMC rep {rep+1} failed: {e}")

        # --- Pathfinder ---
        try:
            t0 = time.time()
            fit_pf = model.pathfinder(
                data=stan_ibd,
                inits=init_dict,
                show_console=False,
            )
            elapsed_pf = time.time() - t0

            all_vars = fit_pf.stan_variables()
            pmeans = {
                name: draws.mean(axis=0)
                for name, draws in all_vars.items()
            }
            results_pf[cm_val].append(pmeans)
            times_pf[cm_val].append(elapsed_pf)

        except Exception as e:
            print(f"    [WARN] Pathfinder rep {rep+1} failed: {e}")


# ====================================================================
# 5. Extract parameters
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
# 6. Plot: 5 subplots (4 bias + 1 time)
# ====================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes_flat = axes.flatten()

param_names = list(true_vals.keys())

for idx, name in enumerate(param_names):
    ax = axes_flat[idx]
    tv = true_vals[name]

    mcmc_bias_means = []
    mcmc_bias_stds = []
    pf_bias_means = []
    pf_bias_stds = []

    for cm_val in cm_values:
        mcmc_vals = np.array([extract_param(pm, name) for pm in results_mcmc[cm_val]])
        pf_vals = np.array([extract_param(pm, name) for pm in results_pf[cm_val]])

        mcmc_bias = (mcmc_vals - tv) / tv * 100
        pf_bias = (pf_vals - tv) / tv * 100

        mcmc_bias_means.append(mcmc_bias.mean())
        mcmc_bias_stds.append(mcmc_bias.std() / np.sqrt(len(mcmc_bias)))
        pf_bias_means.append(pf_bias.mean())
        pf_bias_stds.append(pf_bias.std() / np.sqrt(len(pf_bias)))

    x = np.arange(len(cm_values))
    width = 0.3

    ax.bar(x - width/2, mcmc_bias_means, width, yerr=mcmc_bias_stds,
           color="#378ADD", alpha=0.6, edgecolor="#185FA5", linewidth=1.2,
           label="MCMC", capsize=3)
    ax.bar(x + width/2, pf_bias_means, width, yerr=pf_bias_stds,
           color="#E8A838", alpha=0.6, edgecolor="#C07010", linewidth=1.2,
           label="Pathfinder", capsize=3)

    ax.axhline(0, color="#E24B4A", ls="--", lw=1.5, label="No bias")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{cm}" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_ylabel("Bias (%)")
    ax.set_title(name)
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.15)

# --- Plot 5: Running time ---
ax = axes_flat[4]

mcmc_time_means = []
mcmc_time_stds = []
pf_time_means = []
pf_time_stds = []

for cm_val in cm_values:
    mt = np.array(times_mcmc[cm_val])
    pt = np.array(times_pf[cm_val])
    mcmc_time_means.append(mt.mean())
    mcmc_time_stds.append(mt.std() / np.sqrt(len(mt)))
    pf_time_means.append(pt.mean())
    pf_time_stds.append(pt.std() / np.sqrt(len(pt)))

x = np.arange(len(cm_values))
ax.bar(x - width/2, mcmc_time_means, width, yerr=mcmc_time_stds,
       color="#378ADD", alpha=0.6, edgecolor="#185FA5", linewidth=1.2,
       label="MCMC", capsize=3)
ax.bar(x + width/2, pf_time_means, width, yerr=pf_time_stds,
       color="#E8A838", alpha=0.6, edgecolor="#C07010", linewidth=1.2,
       label="Pathfinder", capsize=3)

ax.set_xticks(x)
ax.set_xticklabels([f"{cm}" for cm in cm_values])
ax.set_xlabel("Genome length (cM)")
ax.set_ylabel("Time (seconds)")
ax.set_title("Running time per replicate")
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.15)

# Hide unused subplot
axes_flat[5].set_visible(False)

plt.suptitle("MCMC vs Pathfinder: Bias and Efficiency", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("bias_efficiency_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: bias_efficiency_comparison.png")

# ====================================================================
# 7. Print summary table
# ====================================================================
print("\n" + "="*80)
print("SUMMARY TABLE")
print("="*80)
print(f"{'cm':>6s} | {'Method':>12s} | {'Ne bias%':>10s} | {'Tadmix bias%':>12s} | "
      f"{'frac bias%':>10s} | {'Tmerge bias%':>12s} | {'Time(s)':>8s}")
print("-"*80)

for cm_val in cm_values:
    for label, res, tm in [
        ("MCMC", results_mcmc, times_mcmc),
        ("Pathfinder", results_pf, times_pf),
    ]:
        biases = {}
        for name, tv in true_vals.items():
            vals = np.array([extract_param(pm, name) for pm in res[cm_val]])
            biases[name] = (vals.mean() - tv) / tv * 100
        mean_time = np.mean(tm[cm_val])
        print(f"{cm_val:6d} | {label:>12s} | {biases['Effective population size']:+9.2f}% | "
              f"{biases['Admixture time']:+11.2f}% | "
              f"{biases['Admixture fraction']:+9.2f}% | "
              f"{biases['Post-admixture merge time']:+11.2f}% | "
              f"{mean_time:8.2f}")
