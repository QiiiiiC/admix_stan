"""
Run 100 replicates of IBD-based demographic inference using:
  - Block-bootstrap resampling (simulate once, resample many times)
  - CmdStan pathfinder (fast approximate posterior)

Produces a violin plot of posterior means across replicates,
one panel per parameter, one violin per genome-length value.
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import msprime
from demography import DemographicTopology
import numpy as np

# ── Your existing imports ──
# from demographic_topology import DemographicTopology  # your topology class
from simulation_methods import simulate_msprime, build_ibd_stan_data
from bootstrap_ibd import (
    calculate_ibd_blocks_mrca,
    resample_ibd_with_bootstrap_variance,
)
from cmdstanpy import CmdStanModel


# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 200
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4          # must match your recombination_rate scaling
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 40, 'b': 40, 'c': 40, 'admix': 40}

# Largest cm value determines how big the simulation genome must be
cm_values = [50, 100, 150, 200, 250, 300, 400, 500]
SIM_CM = max(cm_values)     
SIM_SEQ_LEN = SIM_CM / CM_PER_UNIT   

bins = [
    [0.5, 0.55], [0.55,0.6], [0.6,0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0],[5.0,10.0], [10.0, SIM_CM]
]

# True parameter values (for the plot)
true_vals = {
    "Admixture time": 15,
    "Effective population size": 6000,
    "Admixture fraction": 0.75,
    "Post-admixture merge time": 200,
}

# ====================================================================
# 1. Define topology (same as yours)
# ====================================================================
dem = DemographicTopology(['a', 'admix', 'b', 'c'])
dem.add_admixture_event('admix', 'admixPa1', 'admixPa2')
dem.add_merge_event('a', 'admixPa1', 'aAdmix')
dem.add_merge_event('b', 'admixPa2', 'bAdmix')
dem.add_merge_event('aAdmix', 'bAdmix', 'subAnc')
dem.add_merge_event('subAnc', 'c', 'anc')
for node in ['a','b','c','admix','admixPa1','admixPa2','aAdmix','bAdmix','subAnc','anc']:
    dem.set_node_ne(node, 3000)
dem.set_admixture_parameters('admix', 15, 0.75, 'admixPa1')
dem.set_merge_time('aAdmix', 35)
dem.set_merge_time('bAdmix', 45)
dem.set_merge_time('subAnc', 200)
dem.set_merge_time('anc', 500)
dem.finalize_root()

# ====================================================================
# 2. Simulate ONCE with the largest genome
# ====================================================================
print(f"Simulating {SIM_CM} cM genome...")
mts = simulate_msprime(
    dem,
    sequence_length=SIM_SEQ_LEN,
    recombination_rate=RECOMB_RATE,
    mutation_rate=MUT_RATE,
    samples_per_pop=SAMPLES_PER_POP,
)

# ====================================================================
# 3. Compute block-level IBD (expensive step — done ONCE)
# ====================================================================
print("Computing block-level IBD segments...")
# Note: for the block computation, use bins with the last bin up to SIM_CM
raw_block_ibd, n_blocks, pop_samples, pop_ids = calculate_ibd_blocks_mrca(
    mts,
    bins=bins,
    block_size_cm=BLOCK_SIZE_CM,
    cm_per_unit=CM_PER_UNIT,
)
print(f"Got {n_blocks} blocks of {BLOCK_SIZE_CM} cM each.")

# ====================================================================
# 4. Compile Stan model once
# ====================================================================
model = CmdStanModel(stan_file="ibd_model_Nfixed.stan")

# ====================================================================
# 5. Run replicates
# ====================================================================
# results[cm_val] = list of 100 dicts, each dict has parameter posterior means
results = {cm_val: [] for cm_val in cm_values}

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} replicates)")
    print(f"{'='*60}")

    # Adjust last bin upper bound to match this cm_val
    bins_cm = [list(b) for b in bins]
    bins_cm[-1][1] = cm_val  # last bin goes up to cm_val

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        # --- Block-bootstrap resample ---
        ibd_mean, ibd_var = resample_ibd_with_bootstrap_variance(
            raw_block_ibd,
            n_blocks_total=n_blocks,
            pop_samples=pop_samples,
            pop_ids=pop_ids,
            bins=bins_cm,
            target_cm=cm_val,
            block_size_cm=BLOCK_SIZE_CM,
            n_var_bootstraps=1000,
            rng=rng,
        )

        # --- Build Stan data ---
        stan_ibd = build_ibd_stan_data(
            dem, ibd_mean, ibd_var, bins_cm,
            T_max=100000, cm=cm_val,
        )

        # --- Pathfinder ---
        init_dict = {
            "times": [100.0] * stan_ibd['n_events'],
            "admixture_fractions": [0.5] * stan_ibd['n_admixture'],
            "effective_N": 10000.0,
        }

        try:
            fit = model.pathfinder(
                data=stan_ibd,
                inits=init_dict,
                show_console=False,
            )

            all_vars = fit.stan_variables()
            pmeans = {
                name: draws.mean(axis=0)
                for name, draws in all_vars.items()
            }
            results[cm_val].append(pmeans)

        except Exception as e:
            print(f"    [WARN] replicate {rep+1} failed: {e}")
            # Skip failed fits
            continue

# ====================================================================
# 6. Extract posterior means for plotting
# ====================================================================
def extract_param(pmeans_dict, param_name):
    """Pull a scalar posterior mean from a stan_variables() mean dict."""
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
# 7. Plot
# ====================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, (name, tv) in zip(axes.flatten(), true_vals.items()):
    samples_list = []

    for cm_val in cm_values:
        vals = [
            extract_param(pm, name)
            for pm in results[cm_val]
        ]
        samples_list.append(np.array(vals))

    # Box plot
    bp = ax.boxplot(
        samples_list,
        positions=range(len(cm_values)),
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="#BA7517", linewidth=2),
        whiskerprops=dict(color="k", linewidth=1.2),
        capprops=dict(color="k", linewidth=1.2),
        boxprops=dict(facecolor="#378ADD", alpha=0.3, edgecolor="#185FA5", linewidth=1.2),
    )

    ax.axhline(tv, color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}")
    ax.legend(fontsize=10)
    ax.set_xticks(range(len(cm_values)))
    ax.set_xticklabels([f"{cm} cM" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_title(name)
    ax.grid(axis="y", alpha=0.15)

plt.tight_layout()
plt.savefig("posterior_mean_convergence.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: posterior_mean_convergence.png")
