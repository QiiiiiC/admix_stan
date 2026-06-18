"""Nvarying tree ((a,b),c),d): IBD vs SNP vs Mixed, IBD from the TRUE
tree-sequence MRCA spans (calculate_ibd_blocks_mrca).

Only the parameter relative-error box plot (Mixed model) is produced;
the topology diagram is also saved in this folder.
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")
# Write every figure (results + topology diagram) next to this script.
_HERE = os.path.dirname(os.path.abspath(__file__))

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
from relative_error import plot_relative_error_boxplot
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
N_REPLICATES   = 100
BLOCK_SIZE_CM  = 50.0
CM_PER_UNIT    = 1e-6       # 1 cM per Mb (standard human)
RECOMB_RATE    = 1e-8       # crossovers/bp/generation
MUT_RATE       = 1.25e-8    # mutations/bp/generation
SAMPLES_PER_POP = {'a': 20, 'b': 20, 'c': 20, 'd': 20}
# Haploid sample count per population (diploid samples x ploidy 2)
N_HAPLOID_PER_POP = {_p: 2 * _c for _p, _c in SAMPLES_PER_POP.items()}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [1000]

# Compare bins from 1.0 cM (matches the admix_thru_b_plus_ancient runs)
bins = [
    [1.0, 1.25],[1.25,1.5],
    [1.5, 1.75],[1.75,2.0], [2.0, 2.5],[2.5,3.0],
    [3.0,4.0],[4.0,5.0],[5.0, 6.0],[6.0,7.5],
    [7.5, 12.0],[12.0,20.0], [20.0, BLOCK_SIZE_CM]
]

# Haploid Ne per node (true data-generating values); dem stores diploid Ne (// 2)
NE_HAPLOID = {
    'a': 4000,
    'b': 15000,
    'c': 10000,
    'd': 20000,
    'ab': 10000,
    'abc': 15000,
    'root': 10000,
}

true_vals = {
    "Merge time a+b": 40,
    "Merge time ab+c": 500,
    "Merge time abc+d": 1200,
    "Ne(a)": NE_HAPLOID['a'],
    "Ne(b)": NE_HAPLOID['b'],
    "Ne(c)": NE_HAPLOID['c'],
}

SNP_CUTOFF_TIME = 1200
SNP_MIN_MAF = 0.05

# ====================================================================
# 1. Define topology — pure tree, no admixture
# ====================================================================
dem = DemographicTopology(['a', 'b', 'c', 'd'])
dem.add_merge_event('a', 'b', 'ab')
dem.add_merge_event('ab', 'c', 'abc')
dem.add_merge_event('abc', 'd', 'root')

# Node order: a(0), b(1), c(2), d(3), ab(4), abc(5), root(6)
# msprime uses diploid Ne (haploid // 2)
for node in ['a', 'b', 'c', 'd', 'ab', 'abc', 'root']:
    dem.set_node_ne(node, NE_HAPLOID[node] // 2)

dem.set_merge_time('ab', 20)
dem.set_merge_time('abc', 700)
dem.set_merge_time('root', 1700)
dem.finalize_root()

# --- Save a diagram of the (true) data-generating topology in this folder ---
_topo_fig, _topo_ax = plt.subplots(figsize=(10, 6))
dem.plot_demography(scale=True, ax=_topo_ax)
_topo_ax.set_title('topology_tree (true topology)')
_topo_fig.savefig(os.path.join(_HERE, 'topology_tree.png'), dpi=200, bbox_inches='tight')
plt.close(_topo_fig)
print(f"[topology] saved {os.path.join(_HERE, 'topology_tree.png')}")

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
pop_ids = None

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

# Collect full stan_variables() dicts for the per-model relative-error box plots
# (IBD-only, SNP-only, Mixed; true topology `dem`, largest genome length only).
rel_err_ibd = []
rel_err_snp = []
rel_err_mixed = []
REL_ERR_CM = cm_values[-1]

init_nvarying = {
    "times": [100.0] * n_events,
    "mu_log": float(np.log(15000.0)),
    "sigma_log": 0.3,
    "Ne_raw": [0.0] * n_nodes,
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
                dem, ibd_mean, ibd_var, bins, n_samples_per_pop=N_HAPLOID_PER_POP,
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
            if cm_val == REL_ERR_CM:
                rel_err_ibd.append(all_vars)
        except Exception as e:
            print(f"    [WARN] IBD rep {rep+1} failed: {e}")

        # --- Model 2: SNP-only ---
        try:
            stan_data = build_snp_stan_data(
                dem, w_hat, w_se,
            )
            fit = snp_model.pathfinder(
                data=stan_data, inits=init_nvarying, show_console=False,
                psis_resample=True,
            )
            all_vars = fit.stan_variables()
            pmeans = {name: draws.mean(axis=0) for name, draws in all_vars.items()}
            results_snp[cm_val].append(pmeans)
            elbo_snp[cm_val].append(extract_elbo(fit))
            if cm_val == REL_ERR_CM:
                rel_err_snp.append(all_vars)
        except Exception as e:
            print(f"    [WARN] SNP rep {rep+1} failed: {e}")

        # --- Model 3: Mixed ---
        try:
            stan_data = build_mixed_stan_data(
                dem, ibd_mean, ibd_var, bins, w_hat, w_se, n_samples_per_pop=N_HAPLOID_PER_POP,
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
            if cm_val == REL_ERR_CM:
                rel_err_mixed.append(all_vars)
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
# 6. Print summary table
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

# ====================================================================
# 7. Parameter relative-error box plot (Mixed model, true topology,
#    largest genome length)
# ====================================================================
for _m_label, _m_tag, _m_vars in [
    ("IBD-only", "ibd",   rel_err_ibd),
    ("SNP-only", "snp",   rel_err_snp),
    ("Mixed",    "mixed", rel_err_mixed),
]:
    plot_relative_error_boxplot(
        dem, _m_vars, varying=True,
        outpath=os.path.join(_HERE, f"relative_error_trueMRCA_{_m_tag}.png"),
        title=f"Parameter relative error ({_m_label}, true MRCA IBD)",
    )
