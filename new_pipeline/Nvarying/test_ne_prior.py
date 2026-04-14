"""
Test: how does the Ne gamma prior affect the inferred Ne under T1?

We simulate once from the true pure tree ((a,b),c),d) with Ne(a)=200, Ne(rest)=20000,
then fit the Nvarying SNP model under T1 with several different gamma(shape, rate)
priors on Ne and report the posterior means for Ne(a), Ne(b), Ne(ab), Ne(abc).

True branch-drift ratios (fixed by the data):
    x_a   = 20/200   = 0.1
    x_b   = 20/20000 = 0.001
    x_ab  = 680/20000 = 0.034
    x_abc = 1000/20000 = 0.05
    => Ne(b)/Ne(a) = x_a/x_b = 100
"""

import sys, os
_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt

from simulation_methods import simulate_msprime, build_snp_stan_data
from inference_methods import (
    prepare_snp_blocks,
    pool_snp_blocks,
    resample_snp_covariance,
)
from demography import DemographicTopology
from cmdstanpy import CmdStanModel


# ===== Config =====
CM_PER_UNIT = 1e-4
RECOMB_RATE = 1e-6
MUT_RATE = 1e-6
BLOCK_SIZE_CM = 50.0
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15}
N_SIMS = 20
SIM_CM_EACH = 50
SNP_CUTOFF_TIME = 1700
SNP_MIN_MAF = 0.05

NE_A = 200
NE_REST = 20000


# ===== Build T1 =====
dem = DemographicTopology(['a', 'b', 'c', 'd'])
dem.add_merge_event('a', 'b', 'ab')
dem.add_merge_event('ab', 'c', 'abc')
dem.add_merge_event('abc', 'd', 'root')
dem.set_node_ne('a', NE_A // 2)
for node in ['b', 'c', 'd', 'ab', 'abc', 'root']:
    dem.set_node_ne(node, NE_REST // 2)
dem.set_merge_time('ab', 20)
dem.set_merge_time('abc', 700)
dem.set_merge_time('root', 1700)
dem.finalize_root()


# ===== Simulate =====
print("Simulating...")
snp_blocks_list = []
for sim_i in range(N_SIMS):
    mts = simulate_msprime(
        dem,
        sequence_length=SIM_CM_EACH / CM_PER_UNIT,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )
    snp_blks, _ = prepare_snp_blocks(
        mts, dem,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME,
        min_maf=SNP_MIN_MAF,
    )
    snp_blocks_list.append(snp_blks)

pooled_snp = pool_snp_blocks(snp_blocks_list)
total_blocks = len(pooled_snp)

leaves = dem.initial_leaves
n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])

# Use all pooled blocks for a precise covariance estimate
all_idx = np.arange(total_blocks)
w_hat, w_se = resample_snp_covariance(
    pooled_snp, all_idx,
    n_haploid_per_pop=n_haploid,
    se_block_size=50,
)

print(f"\nw_hat (double-centered) =\n{w_hat}")
print(f"\nw_se =\n{w_se}")


# ===== Figure out Ne indexing =====
node_names = list(dem.nodes.keys())
print(f"\nNode order (Stan index = position + 1): {node_names}")
idx_a   = node_names.index('a')
idx_b   = node_names.index('b')
idx_c   = node_names.index('c')
idx_d   = node_names.index('d')
idx_ab  = node_names.index('ab')
idx_abc = node_names.index('abc')


# ===== Compile model =====
model = CmdStanModel(
    stan_file=os.path.join(_MODELS, "snp_model_Nvarying_priortest.stan")
)

# ===== Prior settings to test =====
# gamma(shape, rate) has mean = shape/rate, mode = (shape-1)/rate (for shape>1)
prior_settings = [
    ("baseline (mean=10000)", 2.0, 0.0002),
    ("lower center (mean=500)", 2.0, 0.004),
    ("higher center (mean=200000)", 2.0, 0.00001),
    ("weak/flat (mean=10000, var large)", 1.0, 0.0001),
    ("tight around 10000 (shape=20)", 20.0, 0.002),
]


# ===== Fit and report =====
n_events_t = len(dem.ordered_events)
n_nodes_t  = len(dem.nodes)

init = {"times": [100.0] * n_events_t, "Ne": [10000.0] * n_nodes_t}

results = []
for label, shape, rate in prior_settings:
    print("\n" + "=" * 70)
    print(f"Prior: {label}  |  gamma(shape={shape}, rate={rate})")
    print(f"  prior mean = {shape/rate:.1f}, mode = {max(0, (shape-1)/rate):.1f}")
    print("=" * 70)

    sd = build_snp_stan_data(dem, w_hat, w_se, effective_N=NE_REST)
    sd['ne_prior_shape'] = shape
    sd['ne_prior_rate']  = rate

    fit = model.pathfinder(
        data=sd, inits=init, show_console=False, psis_resample=True,
        seed=2025,
    )

    all_vars = fit.stan_variables()
    Ne_draws    = all_vars['Ne']       # (ndraws, n_nodes)
    times_draws = all_vars['times']    # (ndraws, n_events)

    ne_mean   = Ne_draws.mean(axis=0)
    ne_median = np.median(Ne_draws, axis=0)
    t_mean    = times_draws.mean(axis=0)

    print(f"  times posterior mean: {t_mean}  (true = [20, 680, 1000])")
    print(f"  Ne(a)   mean={ne_mean[idx_a]:>10.1f}  median={ne_median[idx_a]:>10.1f}  (true 200)")
    print(f"  Ne(b)   mean={ne_mean[idx_b]:>10.1f}  median={ne_median[idx_b]:>10.1f}  (true 20000)")
    print(f"  Ne(c)   mean={ne_mean[idx_c]:>10.1f}  median={ne_median[idx_c]:>10.1f}  (true 20000)")
    print(f"  Ne(d)   mean={ne_mean[idx_d]:>10.1f}  median={ne_median[idx_d]:>10.1f}  (true 20000)")
    print(f"  Ne(ab)  mean={ne_mean[idx_ab]:>10.1f}  median={ne_median[idx_ab]:>10.1f}  (true 20000)")
    print(f"  Ne(abc) mean={ne_mean[idx_abc]:>10.1f}  median={ne_median[idx_abc]:>10.1f}  (true 20000)")
    print(f"  ratio Ne(b)/Ne(a) = {ne_mean[idx_b] / ne_mean[idx_a]:.2f}  (should be ~100)")

    results.append({
        'label': label, 'shape': shape, 'rate': rate,
        'ne_a': ne_mean[idx_a], 'ne_b': ne_mean[idx_b],
        't1': t_mean[0],
        'ne_a_draws': Ne_draws[:, idx_a].copy(),
        'ne_b_draws': Ne_draws[:, idx_b].copy(),
        't1_draws':   times_draws[:, 0].copy(),
    })


# ===== Plot: posterior vs prior =====
n_priors = len(results)
x = np.arange(n_priors)
labels_short = [f"mean={r['shape']/r['rate']:.0f}\nshape={r['shape']}"
                for r in results]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

def violin_panel(ax, draws_list, truth, title, ylabel, log=False):
    vp = ax.violinplot(
        draws_list, positions=x, widths=0.7,
        showmeans=True, showmedians=True, showextrema=False,
    )
    for body in vp['bodies']:
        body.set_facecolor("#4CAF50")
        body.set_alpha(0.5)
        body.set_edgecolor("k")
        body.set_linewidth(0.5)
    vp['cmeans'].set_color("k")
    vp['cmedians'].set_color("k")
    vp['cmedians'].set_linestyle("--")
    if truth is not None:
        ax.axhline(truth, color='red', ls='--', lw=1.2,
                   label=f"truth = {truth:g}")
        ax.legend(fontsize=9, loc='best')
    ax.set_xticks(x)
    ax.set_xticklabels(labels_short, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", alpha=0.15)
    if log:
        ax.set_yscale("log")

# Ne(a)
violin_panel(
    axes[0, 0],
    [r['ne_a_draws'] for r in results],
    truth=NE_A,
    title="Ne(a) posterior across priors",
    ylabel="Ne(a)",
    log=True,
)

# Ne(b)
violin_panel(
    axes[0, 1],
    [r['ne_b_draws'] for r in results],
    truth=NE_REST,
    title="Ne(b) posterior across priors",
    ylabel="Ne(b)",
    log=True,
)

# t_1
violin_panel(
    axes[1, 0],
    [r['t1_draws'] for r in results],
    truth=20.0,
    title="times[1] = t_ab posterior across priors",
    ylabel="t_1",
    log=True,
)

# ratio
ratio_draws = [r['ne_b_draws'] / r['ne_a_draws'] for r in results]
violin_panel(
    axes[1, 1],
    ratio_draws,
    truth=100.0,
    title="Ne(b) / Ne(a) ratio across priors",
    ylabel="ratio",
    log=True,
)

fig.suptitle(
    f"Nvarying SNP model: Ne inference under T1 vs. Ne gamma prior\n"
    f"(pooled genome = {N_SIMS * SIM_CM_EACH} cM, true Ne(a)={NE_A}, Ne(rest)={NE_REST})",
    fontsize=13, fontweight="bold"
)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("test_ne_prior.png", dpi=180, bbox_inches="tight")
print("\nSaved: test_ne_prior.png")


# ===== Summary =====
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"{'prior':<38s} {'t1':>8s} {'Ne(a)':>10s} {'Ne(b)':>10s} {'ratio':>8s}")
print("-" * 80)
for r in results:
    print(f"{r['label']:<38s} {r['t1']:>8.1f} {r['ne_a']:>10.1f} {r['ne_b']:>10.1f} "
          f"{r['ne_b']/r['ne_a']:>8.1f}")
print(f"{'TRUTH':<38s} {20.0:>8.1f} {200.0:>10.1f} {20000.0:>10.1f} {100.0:>8.1f}")
