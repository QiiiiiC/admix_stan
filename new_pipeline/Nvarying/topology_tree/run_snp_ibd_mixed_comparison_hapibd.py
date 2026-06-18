"""Nvarying tree ((a,b),c),d): IBD vs SNP vs Mixed, with IBD segments from
Hap-IBD (run on a VCF + uniform PLINK map) instead of true MRCA spans.
Mirrors the Nfixed Hap-IBD variant. Only the parameter relative-error box
plot (Mixed model) is produced; the topology diagram is also saved here.

Requires hap-ibd in PATH (or HAPIBD_CMD) and Java >= 8.
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
import gzip, shutil, subprocess, time, glob
from pathlib import Path

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
# 0b. Hap-IBD helpers (IBD detection) -- verbatim from the Nfixed variant
# ====================================================================
# Hap-IBD knobs --- the segment-length floor matches the smallest IBD bin (1.0 cM).
HAPIBD_MIN_SEED    = 1.0
HAPIBD_MIN_OUTPUT  = 1.0
HAPIBD_THREADS     = int(os.environ.get("HAPIBD_THREADS", "1"))


def get_hapibd_cmd_template():
    env = os.environ.get("HAPIBD_CMD")
    if env:
        return env
    if shutil.which("hap-ibd") is None:
        raise RuntimeError(
            "hap-ibd not found in PATH. "
            "Set HAPIBD_CMD, e.g.:\n"
            "  export HAPIBD_CMD='hap-ibd gt={vcf_path} map={map_path} "
            "out={output_prefix} min-seed={min_seed} min-output={min_output} "
            "nthreads={threads}'"
        )
    return (
        "hap-ibd gt={vcf_path} map={map_path} out={output_prefix} "
        "min-seed={min_seed} min-output={min_output} nthreads={threads}"
    )

def run_cmd(cmd, quiet=True):
    if not quiet:
        print(f"  Running: {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print("STDOUT:", res.stdout)
        print("STDERR:", res.stderr)
        raise RuntimeError(f"Command failed: {cmd}")

def write_vcf_and_map(ts, output_dir: Path, prefix: str):
    """Write VCF directly from msprime + a 2-line uniform PLINK map.

    Assumes msprime was run with sim_ancestry(..., discrete_genome=True),
    which makes all mutation positions integer bp from the start. No
    position transform needed.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    vcf_path = output_dir / f"{prefix}.vcf"
    map_path = output_dir / f"{prefix}.map"

    with open(vcf_path, "w") as f:
        ts.write_vcf(f, contig_id="1")

    seq_len_bp = int(ts.sequence_length)
    total_cm   = seq_len_bp * CM_PER_UNIT
    with open(map_path, "w") as f:
        f.write(f"1\t.\t0.000000\t1\n")
        f.write(f"1\t.\t{total_cm:.6f}\t{seq_len_bp}\n")

    # Read back sample names from the VCF header
    sample_names = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                sample_names = line.strip().split()[9:]
                break
    if not sample_names:
        raise ValueError(f"No sample columns in VCF header: {vcf_path}")

    return vcf_path, map_path, sample_names

def find_hapibd_output(output_prefix: Path):
    """Hap-IBD writes <prefix>.ibd.gz by default."""
    for ext in [".ibd.gz", ".ibd"]:
        p = Path(f"{output_prefix}{ext}")
        if p.exists():
            return p
    candidates = list(output_prefix.parent.glob(output_prefix.name + "*"))
    for p in candidates:
        if p.name.endswith(".ibd.gz") or p.suffix == ".ibd":
            return p
    raise FileNotFoundError(
        f"Hap-IBD output not found for prefix {output_prefix}. "
        f"Files in dir: {list(output_prefix.parent.iterdir())}"
    )

def parse_hapibd_segments(output_file: Path, sample_to_node):
    """Read Hap-IBD .ibd / .ibd.gz output.
    Format per line: sample1 hap1 sample2 hap2 chrom start_bp end_bp length_cm
    We map (sample, hap) -> tskit haplotype node ID.
    """
    opener = gzip.open if str(output_file).endswith(".gz") else open
    segments = []  # list of (node_u, node_v, length_cm)
    with opener(output_file, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.strip().split()
            if len(toks) < 8:
                continue
            s1, h1, s2, h2 = toks[0], toks[1], toks[2], toks[3]
            try:
                length_cm = float(toks[7])
            except ValueError:
                continue

            key1 = (s1, int(h1))
            key2 = (s2, int(h2))
            if key1 not in sample_to_node or key2 not in sample_to_node:
                continue
            segments.append((sample_to_node[key1], sample_to_node[key2], length_cm))
    return segments

def build_sample_to_node(ts, sample_names):
    """Map (vcf_sample_name, vcf_hap_index_1_or_2) -> tskit sample node ID.
    msprime's write_vcf emits diploid samples (tsk_0, tsk_1, ...) and
    pairs consecutive sample nodes into the two haplotypes of each diploid.
    """
    sample_nodes = list(ts.samples())
    if len(sample_nodes) != 2 * len(sample_names):
        raise ValueError(
            f"Sample-node count mismatch: {len(sample_nodes)} nodes vs "
            f"{len(sample_names)} diploid samples in VCF"
        )
    mapping = {}
    for i, name in enumerate(sample_names):
        mapping[(name, 1)] = sample_nodes[2 * i]
        mapping[(name, 2)] = sample_nodes[2 * i + 1]
    return mapping

def hapibd_segments_to_packed(ts, segments, bins, block_size_cm, cm_per_unit):
    """Convert a list of (node_u, node_v, length_cm) into the
    packed[(b_i, i, j)] = (num_pairs, n_blocks) format.

    For one simulation of length SIM_CM_EACH = block_size_cm, n_blocks = 1
    and each segment contributes its full length / block_size_cm (the
    'fraction of the block covered by the segment') to that single block.

    Returns
    -------
    packed : dict   keyed by (b_i, i, j)
    n_blocks : int  (always 1 for the matched config)
    pop_samples : dict {pop_id: [sample_node_ids]}
    pop_ids : np.ndarray
    pair_info : dict  (same schema as calculate_ibd_blocks_mrca)
    """
    sample_nodes = ts.samples()
    node_to_pop = ts.nodes_population[sample_nodes]
    pop_ids = np.unique(node_to_pop)
    num_pops = len(pop_ids)

    pop_samples = defaultdict(list)
    for u in sample_nodes:
        pop_samples[node_to_pop[u]].append(u)

    genome_length_cm = ts.sequence_length * cm_per_unit
    n_blocks = int(np.floor(genome_length_cm / block_size_cm))
    if n_blocks == 0:
        raise ValueError(
            f"Genome ({genome_length_cm:.2f} cM) shorter than block ({block_size_cm} cM)"
        )

    # --- enumerate pairs in the same deterministic order as the true-IBD code ---
    pair_info = {}
    pair_key_to_index = {}     # (b_i, i, j) -> {(min_node, max_node): row_idx}
    packed_shape = {}          # (b_i, i, j) -> num_pairs

    for i in range(num_pops):
        for j in range(i, num_pops):
            haps_i = pop_samples[pop_ids[i]]
            haps_j = pop_samples[pop_ids[j]]
            n_i, n_j = len(haps_i), len(haps_j)

            pair_a, pair_b, pair_keys = [], [], []
            if i == j:
                for a in range(n_i):
                    for b in range(a + 1, n_i):
                        pair_a.append(a); pair_b.append(b)
                        pair_keys.append(tuple(sorted((haps_i[a], haps_i[b]))))
            else:
                for a in range(n_i):
                    for b in range(n_j):
                        pair_a.append(a); pair_b.append(b)
                        pair_keys.append(tuple(sorted((haps_i[a], haps_j[b]))))

            pair_info[(i, j)] = {
                'pair_a':   np.array(pair_a, dtype=np.int32),
                'pair_b':   np.array(pair_b, dtype=np.int32),
                'n_haps_i': n_i,
                'n_haps_j': n_j,
                'is_within': i == j,
            }

            idx_map = {pk: row for row, pk in enumerate(pair_keys)}
            for b_i_idx in range(len(bins)):
                key = (b_i_idx, i, j)
                pair_key_to_index[key] = idx_map
                packed_shape[key] = len(pair_keys)

    # --- allocate packed arrays ---
    packed = {key: np.zeros((packed_shape[key], n_blocks))
              for key in packed_shape}

    # --- which (i, j) does a given (u, v) belong to? ---
    node_to_pop_idx = {sample_nodes[k]: np.searchsorted(pop_ids, node_to_pop[k])
                       for k in range(len(sample_nodes))}

    # --- assign each Hap-IBD segment to its bin and pair-row ---
    n_segs_kept = 0
    n_segs_skipped = 0
    for u, v, length_cm in segments:
        # bin lookup (matches calculate_ibd_blocks_mrca: min_len < L <= max_len)
        target_bin = -1
        for b_i_idx, (lo, hi) in enumerate(bins):
            if lo < length_cm <= hi:
                target_bin = b_i_idx
                break
        if target_bin == -1:
            n_segs_skipped += 1
            continue

        pi = node_to_pop_idx[u]
        pj = node_to_pop_idx[v]
        i, j = (pi, pj) if pi <= pj else (pj, pi)

        pk = tuple(sorted((u, v)))
        idx_map = pair_key_to_index[(target_bin, i, j)]
        if pk not in idx_map:
            n_segs_skipped += 1
            continue
        row = idx_map[pk]

        # n_blocks = 1 by construction in this experiment; segment contributes
        # length_cm / block_size_cm (fraction of the block) to that one block.
        # If you ever simulate longer than block_size_cm, you'd need to split
        # the segment across block boundaries here.
        frac = length_cm / block_size_cm
        packed[(target_bin, i, j)][row, 0] += frac
        n_segs_kept += 1

    return packed, n_blocks, dict(pop_samples), pop_ids, pair_info, n_segs_kept, n_segs_skipped


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

hapibd_cmd_template = get_hapibd_cmd_template()
print(f"Hap-IBD template: {hapibd_cmd_template}")
work_root = Path(_PARENT) / "hapibd_work"
work_root.mkdir(exist_ok=True)
# Per-script tag: isolates each script's per-simulation dirs under the
# shared work root, so two different run scripts can run concurrently
# without clobbering each other's sim_* inputs/outputs.
_RUN_TAG = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]

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

    sim_dir = work_root / f"{_RUN_TAG}_sim_{sim_i}"
    vcf_path, map_path, sample_names = write_vcf_and_map(mts, sim_dir, "sim")
    hapibd_prefix = sim_dir / "hapibd_out"
    cmd = hapibd_cmd_template.format(
        vcf_path=str(vcf_path), map_path=str(map_path),
        output_prefix=str(hapibd_prefix),
        min_seed=HAPIBD_MIN_SEED, min_output=HAPIBD_MIN_OUTPUT, threads=HAPIBD_THREADS,
    )
    run_cmd(cmd, quiet=True)
    hapibd_out = find_hapibd_output(hapibd_prefix)
    sample_to_node = build_sample_to_node(mts, sample_names)
    segments = parse_hapibd_segments(hapibd_out, sample_to_node)
    packed, n_blocks, pop_samples, pids, pi, n_kept, n_skip = hapibd_segments_to_packed(
        mts, segments, bins, block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
    )
    print(f"  [Hap-IBD] segs kept={n_kept} skipped={n_skip}")
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi
        pop_ids = pids

    snp_blks, _ = prepare_snp_blocks(
        mts, dem,
        block_size_cm=BLOCK_SIZE_CM,
        cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME,
        min_maf=SNP_MIN_MAF,
    )
    snp_blocks_list.append(snp_blks)

    shutil.rmtree(sim_dir, ignore_errors=True)

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
        outpath=os.path.join(_HERE, f"relative_error_hapibd_{_m_tag}.png"),
        title=f"Parameter relative error ({_m_label}, Hap-IBD)",
    )
