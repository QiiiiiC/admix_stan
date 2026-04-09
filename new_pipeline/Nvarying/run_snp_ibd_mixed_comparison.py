"""
Compare three Nvarying models on a pure tree (no admixture):
  1. IBD-only   (ibd_model_Nvarying.stan, per-node Ne)
  2. SNP-only   (snp_model_Nvarying.stan, per-node Ne)
  3. Mixed      (mixed_model_Nvarying.stan, per-node Ne)

Topology: ((a,b),c),d)
  a + b → ab   at t = 20
  ab + c → abc at t = 700
  abc + d → root at t = 1700

Ne(a) = 200 (haploid), all other nodes = 20000 (haploid).
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

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

# ====================================================================
# 0. Configuration
# ====================================================================
N_REPLICATES = 200
BLOCK_SIZE_CM = 50.0
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 250, 300, 400, 500]

bins = [
    [0.5, 0.55], [0.55, 0.6], [0.6, 0.65],[0.65, 0.7], [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    [1.0, 1.5], [1.5, 2.0], [2.0, 5.0], [5.0, 8.0],[8.0,20.0], [20.0, BLOCK_SIZE_CM]
]

NE_A = 200       # haploid Ne for node a
NE_REST = 20000  # haploid Ne for all other nodes

true_vals = {
    "Merge time a+b": 20,
    "Merge time ab+c": 700,
    "Merge time abc+d": 1700,
    "Ne(a)": NE_A,
    "Ne(b)": NE_REST,
    "Ne(c)": NE_REST,
}

SNP_CUTOFF_TIME = 1700
SNP_MIN_MAF = 0.05
FIXED_NE = NE_REST  # passed to build_snp_stan_data (ignored by Nvarying model)

# ====================================================================
# 1. Define topology — pure tree, no admixture
# ====================================================================
dem = DemographicTopology(['a', 'b', 'c', 'd'])
dem.add_merge_event('a', 'b', 'ab')
dem.add_merge_event('ab', 'c', 'abc')
dem.add_merge_event('abc', 'd', 'root')

# Node order: a(0), b(1), c(2), d(3), ab(4), abc(5), root(6)
dem.set_node_ne('a', NE_A // 2)       # msprime uses diploid Ne
for node in ['b', 'c', 'd', 'ab', 'abc', 'root']:
    dem.set_node_ne(node, NE_REST // 2)

dem.set_merge_time('ab', 20)
dem.set_merge_time('abc', 700)
dem.set_merge_time('root', 1700)
dem.finalize_root()

# Node order: a(0), b(1), c(2), d(3), ab(4), abc(5), root(6)
NODE_A_IDX = 0
NODE_B_IDX = 1
NODE_C_IDX = 2

n_nodes = len(dem.nodes)
n_events = len(dem.ordered_events)
n_admix = dem.n_admix  # 0

# ====================================================================
# 2. Simulate and pool blocks
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

    snp_blks, _ = prepare_snp_blocks(
        mts, dem,
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
ibd_model = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nvarying.stan"))
snp_model = CmdStanModel(stan_file=os.path.join(_MODELS, "snp_model_Nvarying.stan"))
mixed_model = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nvarying.stan"))

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

init_nvarying = {
    "times": [100.0] * n_events,
    "Ne": [10000.0] * n_nodes,
}
init_mixed_nvarying = {
    **init_nvarying,
    "kappa_snp": 1.0,
}

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  ({N_REPLICATES} replicates)")
    print(f"{'='*60}")

    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

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

        # --- Resample SNP ---
        leaves = dem.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid,
            se_block_size=50,
        )

        # --- Model 1: IBD-only ---
        try:
            stan_data = build_ibd_stan_data(
                dem, ibd_mean, ibd_var, bins,
                T_max=100000, cm=cm_val,
            )
            fit = ibd_model.pathfinder(
                data=stan_data, inits=init_nvarying, show_console=False,
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
                data=stan_data, inits=init_nvarying, show_console=False,
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
                data=stan_data, inits=init_mixed_nvarying, show_console=False,
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
def extract_param(pmeans_dict, param_name):
    t = pmeans_dict['times']
    ne = pmeans_dict['Ne']
    if param_name == "Merge time a+b":
        return t[0] if np.ndim(t) > 0 else float(t)
    elif param_name == "Merge time ab+c":
        return float(np.sum(t[:2]))
    elif param_name == "Merge time abc+d":
        return float(np.sum(t[:3]))
    elif param_name == "Ne(a)":
        return float(ne[NODE_A_IDX])
    elif param_name == "Ne(b)":
        return float(ne[NODE_B_IDX])
    elif param_name == "Ne(c)":
        return float(ne[NODE_C_IDX])
    else:
        raise ValueError(f"Unknown param: {param_name}")


# ====================================================================
# 6a. Plot 1: Parameter comparison (2x2)
# ====================================================================
fig1, axes1 = plt.subplots(2, 3, figsize=(22, 12))

model_labels = ["IBD-only", "SNP-only", "Mixed"]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]
model_results = [results_ibd, results_snp, results_mixed]
from matplotlib.patches import Patch

for ax, (name, tv) in zip(axes1.flatten(), true_vals.items()):
    n_models = len(model_labels)
    n_cms = len(cm_values)

    positions_all = []
    data_all = []
    colors_all = []

    for cm_idx, cm_val in enumerate(cm_values):
        group_center = cm_idx * (n_models + 1)
        for m_idx in range(n_models):
            pos = group_center + m_idx
            positions_all.append(pos)
            colors_all.append(model_colors[m_idx])

            vals = [extract_param(pm, name) for pm in model_results[m_idx][cm_val]]
            data_all.append(np.array(vals))

    vp = ax.violinplot(
        data_all,
        positions=positions_all,
        widths=0.6,
        showmeans=True,
        showmedians=True,
        showextrema=False,
    )

    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors_all[i])
        body.set_alpha(0.5)
        body.set_edgecolor("k")
        body.set_linewidth(0.5)
    vp['cmeans'].set_color("k")
    vp['cmeans'].set_linewidth(1)
    vp['cmedians'].set_color("k")
    vp['cmedians'].set_linewidth(1)
    vp['cmedians'].set_linestyle("--")

    ax.axhline(tv, color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}")

    group_centers = [cm_idx * (n_models + 1) + (n_models - 1) / 2
                     for cm_idx in range(n_cms)]
    ax.set_xticks(group_centers)
    ax.set_xticklabels([f"{cm}" for cm in cm_values])
    ax.set_xlabel("Genome length (cM)")
    ax.set_title(name)
    ax.grid(axis="y", alpha=0.15)

    handles = [Patch(facecolor=c, alpha=0.5, label=l)
               for c, l in zip(model_colors, model_labels)]
    handles.append(plt.Line2D([0], [0], color="#E24B4A", ls="--", lw=1.5, label=f"True = {tv}"))
    ax.legend(handles=handles, fontsize=8)

fig1.suptitle("Nvarying (no admixture): IBD vs SNP vs Mixed\n"
              f"Tree: ((a,b),c),d)  |  Ne(a)={NE_A}, Ne(rest)={NE_REST}\n"
              f"Merge times: a+b={20}, ab+c={700}, abc+d={1700}",
              fontsize=12, fontweight="bold")
fig1.tight_layout()
fig1.savefig("snp_ibd_mixed_comparison.png", dpi=200, bbox_inches="tight")
plt.show()
print("Saved: snp_ibd_mixed_comparison.png")

# ====================================================================
# 6b. Print summary table
# ====================================================================
print("\n" + "=" * 100)
print("SUMMARY TABLE — Bias (%) of posterior mean  [Nvarying, no admixture]")
print("=" * 100)
print(f"{'cm':>6s} | {'Model':>10s} | {'T(a+b)%':>9s} | {'T(ab+c)%':>10s} | "
      f"{'T(abc+d)%':>11s} | {'Ne(a)%':>9s} | {'Ne(b)%':>9s} | {'Ne(c)%':>9s} | {'n_ok':>5s}")
print("-" * 120)

for cm_val in cm_values:
    for label, res in [
        ("IBD", results_ibd),
        ("SNP", results_snp),
        ("Mixed", results_mixed),
    ]:
        biases = {}
        for pname, tv in true_vals.items():
            vals = np.array([extract_param(pm, pname) for pm in res[cm_val]])
            biases[pname] = (np.nanmean(vals) - tv) / tv * 100

        print(f"{cm_val:6d} | {label:>10s} | "
              f"{biases['Merge time a+b']:+8.2f}% | "
              f"{biases['Merge time ab+c']:+9.2f}% | "
              f"{biases['Merge time abc+d']:+10.2f}% | "
              f"{biases['Ne(a)']:+8.2f}% | "
              f"{biases['Ne(b)']:+8.2f}% | "
              f"{biases['Ne(c)']:+8.2f}% | "
              f"{len(res[cm_val]):5d}")
