"""
Hap-IBD version of the Nfixed admix-through-b + ancient-admix experiment.

Mirrors run_admix_thru_b_plus_ancient.py exactly, except IBD segments come
from Hap-IBD (run on a VCF + uniform PLINK map) instead of from the true
tree-sequence MRCA spans.

Pipeline (per simulation, identical structure to the true-IBD version):
  1. msprime simulate 50 cM at human-realistic params -> tree sequence
  2. Write VCF + 2-line uniform PLINK map (1 cM per Mb)
  3. Run Hap-IBD
  4. Parse segments, pack into packed[(b_i, i, j)] of shape (num_pairs, 1)
     in the same deterministic pair order as calculate_ibd_blocks_mrca.
  5. SNP blocks computed exactly as before.

After all sims, pool with pool_multiple_simulations and resample with
resample_ibd_with_jackknife_variance --- both unchanged.

Requirements:
  - hap-ibd available in PATH (or set HAPIBD_CMD env var)
  - Java >= 8 (Hap-IBD is a JAR)
"""

import sys, os, time, gzip, shutil, subprocess, re
from pathlib import Path
from collections import defaultdict

_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
    pool_multiple_simulations,
    resample_ibd_with_jackknife_variance,
)
from demography import DemographicTopology
from cmdstanpy import CmdStanModel


# ====================================================================
# 0. Configuration  (human-realistic units so Hap-IBD just works)
# ====================================================================
N_REPLICATES   = 100
BLOCK_SIZE_CM  = 50.0
CM_PER_UNIT    = 1e-6       # 1 cM per Mb (standard human)
RECOMB_RATE    = 1e-8       # crossovers/bp/generation
MUT_RATE       = 1.25e-8    # mutations/bp/generation
SAMPLES_PER_POP = {'a': 15, 'b': 15, 'c': 15, 'd': 15}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 300, 500, 750, 1000]


bins = [
    # [0.5, 0.55], [0.55, 0.6], [0.6, 0.65], [0.65, 0.7],
    # [0.7, 0.8], [0.8, 0.9], [0.9, 1.0],
    # [1.0, 1.25],[1.25,1.5], 
    [1.5, 1.75],[1.75,2.0], [2.0, 2.5],[2.5,3.0],
    [3.0,4.0],[4.0,5.0],[5.0, 6.0],[6.0,7.5],
    [7.5, 12.0],[12.0,20.0], [20.0, BLOCK_SIZE_CM]
]

# Single uniform Ne everywhere (haploid; demography API takes diploid = haploid//2)
NE_UNIFORM = 10000

# --- Recent admix on a (through b) ---
T_LOOP_SPLIT  = 10
T_APB_MERGE   = 55
T_LOOP_CLOSE  = 75
ALPHA_LOOP    = 0.5

# --- Ancient admix on c ---
T_OLD = 500
F_OLD = 0.70

# --- Other event times ---
T_LEFT_MERGE  = 700
T_RIGHT_MERGE = 900
T_ROOT        = 1000

SNP_CUTOFF_TIME = T_ROOT
SNP_MIN_MAF = 0.05
FIXED_NE = NE_UNIFORM

# Hap-IBD knobs --- min-output is the segment-length floor in cM
HAPIBD_MIN_SEED    = 1.0     # internal seed threshold (cM)
HAPIBD_MIN_OUTPUT  = 1.0      # report all segments >= this cM
HAPIBD_THREADS     = int(os.environ.get("HAPIBD_THREADS", "1"))


# ====================================================================
# 1. Hap-IBD helpers
# ====================================================================
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


# ====================================================================
# 2. Pack Hap-IBD segments into the same format as calculate_ibd_blocks_mrca
# ====================================================================
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
# 3. Stan helpers (unchanged from true-IBD version)
# ====================================================================
def extract_elbo(fit):
    best = -np.inf
    for stdout_file in fit._runset.stdout_files:
        try:
            with open(stdout_file, 'r') as f:
                for line in f:
                    m = re.search(r'Best Iter:.*?ELBO \(([-+\d.eE]+)\)', line)
                    if m:
                        v = float(m.group(1))
                        if v > best:
                            best = v
        except (OSError, ValueError):
            continue
    return best if np.isfinite(best) else np.nan


def get_admix_event_mapping(dem):
    mapping = []
    admix_counter = 0
    for i, ev in enumerate(dem.ordered_events):
        if ev['type'] == 'ADMIXTURE':
            mapping.append((ev['child'], i, admix_counter))
            admix_counter += 1
    return mapping


def get_merge_event_index(dem, parent_name):
    for i, ev in enumerate(dem.ordered_events):
        if ev['type'] == 'MERGE' and ev['parent'] == parent_name:
            return i
    return None


# ====================================================================
# 4. Topology builders (verbatim from the true-IBD script)
# ====================================================================
def _set_uniform_ne(dem, ne_haploid):
    for n in dem.nodes:
        dem.set_node_ne(n, ne_haploid // 2)


def build_t_true():
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP2', 'b', 'bP')
    d.add_merge_event('bP', 'aP1', 'ab')
    d.add_admixture_event('c', 'cP1', 'cP2')
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d', 'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.set_admixture_parameters('a', time=T_LOOP_SPLIT,
                               fraction_parent_1=ALPHA_LOOP, parent_1_name='aP1')
    d.set_merge_time('bP', T_APB_MERGE)
    d.set_merge_time('ab', T_LOOP_CLOSE)
    d.set_admixture_parameters('c', time=T_OLD,
                               fraction_parent_1=F_OLD, parent_1_name='cP1')
    d.set_merge_time('left', T_LEFT_MERGE)
    d.set_merge_time('right', T_RIGHT_MERGE)
    d.set_merge_time('root', T_ROOT)
    d.finalize_root()
    return d


def build_t_alt1():
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_merge_event('a', 'b', 'ab')
    d.add_admixture_event('c', 'cP1', 'cP2')
    d.add_merge_event('ab', 'cP1', 'left')
    d.add_merge_event('d', 'cP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.finalize_root()
    return d


def build_t_alt2():
    d = DemographicTopology(['a', 'b', 'c', 'd'])
    d.add_admixture_event('a', 'aP1', 'aP2')
    d.add_merge_event('aP2', 'b', 'bP')
    d.add_merge_event('bP', 'aP1', 'ab')
    d.add_merge_event('ab', 'd', 'abd')
    d.add_merge_event('abd', 'c', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.finalize_root()
    return d


# ====================================================================
# 5. Main pipeline
# ====================================================================
print("=" * 60)
print("Hap-IBD estimation pipeline -- Nfixed admix-through-b + ancient")
print(f"  Sim length: {SIM_CM_EACH} cM  ({SIM_SEQ_LEN} bp)")
print(f"  N sims    : {N_SIMS}")
print(f"  CM per bp : {CM_PER_UNIT}  (recomb {RECOMB_RATE}/bp/gen)")
print("=" * 60)

dem_true = build_t_true()
dem_alt1 = build_t_alt1()
dem_alt2 = build_t_alt2()
topology_dems = [dem_true, dem_alt1, dem_alt2]
topology_labels = [
    "T_true (recent thru b + ancient)",
    "T_alt1 (ancient only)",
    "T_alt2 (recent only)",
]
for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, "
          f"n_admix={d.n_admix}, n_nodes={len(d.nodes)}")
n_topos = len(topology_dems)
dem_sim = dem_true

hapibd_cmd_template = get_hapibd_cmd_template()
print(f"\nHap-IBD template: {hapibd_cmd_template}")


# ---- Per-simulation: msprime -> VCF -> Hap-IBD -> packed ----
packed_list     = []
n_blocks_list   = []
snp_blocks_list = []
pair_info       = None
pop_ids         = None

work_root = Path(_PARENT) / "hapibd_work"
work_root.mkdir(exist_ok=True)

for sim_i in range(N_SIMS):
    print(f"\n--- Simulation {sim_i + 1}/{N_SIMS} ({SIM_CM_EACH} cM) ---")
    t0 = time.time()

    mts = simulate_msprime(
        dem_sim,
        sequence_length=SIM_SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        mutation_rate=MUT_RATE,
        samples_per_pop=SAMPLES_PER_POP,
        seed=42 + sim_i,
    )
    t1 = time.time()
    print(f"  [msprime] {t1 - t0:.2f}s | n_trees={mts.num_trees} n_sites={mts.num_sites}")

    sim_dir = work_root / f"sim_{sim_i}"
    vcf_path, map_path, sample_names = write_vcf_and_map(mts, sim_dir, "sim")

    hapibd_prefix = sim_dir / "hapibd_out"
    cmd = hapibd_cmd_template.format(
        vcf_path=str(vcf_path),
        map_path=str(map_path),
        output_prefix=str(hapibd_prefix),
        min_seed=HAPIBD_MIN_SEED,
        min_output=HAPIBD_MIN_OUTPUT,
        threads=HAPIBD_THREADS,
    )
    run_cmd(cmd, quiet=True)
    hapibd_out = find_hapibd_output(hapibd_prefix)
    t2 = time.time()
    print(f"  [Hap-IBD] {t2 - t1:.2f}s -> {hapibd_out.name}")

    sample_to_node = build_sample_to_node(mts, sample_names)
    segments = parse_hapibd_segments(hapibd_out, sample_to_node)

    packed, n_blocks, pop_samples, pids, pi, n_kept, n_skip = \
        hapibd_segments_to_packed(
            mts, segments, bins,
            block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
        )
    t3 = time.time()
    print(f"  [Pack ] {t3 - t2:.2f}s | segs kept={n_kept} skipped={n_skip} "
          f"(out-of-bin or unknown pair)")
    packed_list.append(packed)
    n_blocks_list.append(n_blocks)
    if pair_info is None:
        pair_info = pi
        pop_ids = pids

    # SNP blocks: exactly as in the true-IBD version
    snp_blks, _ = prepare_snp_blocks(
        mts, dem_sim,
        block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
        cutoff_time=SNP_CUTOFF_TIME, min_maf=SNP_MIN_MAF,
    )
    t4 = time.time()
    print(f"  [SNP  ] {t4 - t3:.2f}s")
    snp_blocks_list.append(snp_blks)

    # Free disk
    shutil.rmtree(sim_dir, ignore_errors=True)

pooled_ibd, total_blocks = pool_multiple_simulations(packed_list, n_blocks_list)
pooled_snp = pool_snp_blocks(snp_blocks_list)
print(f"\nTotal pool: {total_blocks} blocks of {BLOCK_SIZE_CM} cM = "
      f"{total_blocks * BLOCK_SIZE_CM} cM from {N_SIMS} simulations")


# ====================================================================
# 6. Stan models (Nfixed)
# ====================================================================
ibd_stan   = CmdStanModel(stan_file=os.path.join(_MODELS, "ibd_model_Nfixed.stan"))
snp_stan   = CmdStanModel(stan_file=os.path.join(_MODELS, "snp_model_Nfixed.stan"))
mixed_stan = CmdStanModel(stan_file=os.path.join(_MODELS, "mixed_model_Nfixed.stan"))

model_labels = ["IBD-only", "SNP-only", "Mixed"]
model_stans  = [ibd_stan, snp_stan, mixed_stan]
model_colors = ["#378ADD", "#4CAF50", "#E8A838"]


# ====================================================================
# 7. Per-topology event indices
# ====================================================================
admix_mappings = [get_admix_event_mapping(d) for d in topology_dems]
bP_close_idx   = [get_merge_event_index(d, 'bP') for d in topology_dems]
ab_close_idx   = [get_merge_event_index(d, 'ab') for d in topology_dems]
for ti, (lbl, mp, bpi, abi) in enumerate(zip(
        topology_labels, admix_mappings, bP_close_idx, ab_close_idx)):
    print(f"  {lbl} admix mapping: {mp}, bP idx: {bpi}, ab idx: {abi}")


# ====================================================================
# 8. Replicates: resample blocks -> fit Stan for each (cm, rep, topo, model)
# ====================================================================
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}
admix_params = {
    label: {ti: {cm: [] for cm in cm_values} for ti in range(n_topos)}
    for label in model_labels
}

rng = np.random.default_rng(seed=2025)

for cm_val in cm_values:
    print(f"\n{'='*60}")
    print(f"  cm = {cm_val} cM  "
          f"({N_REPLICATES} reps x {n_topos} topos x {len(model_labels)} models)")
    print(f"{'='*60}")
    n_blocks_needed = max(1, int(round(cm_val / BLOCK_SIZE_CM)))

    for rep in range(N_REPLICATES):
        if (rep + 1) % 10 == 0:
            print(f"  replicate {rep + 1}/{N_REPLICATES}")

        replace = n_blocks_needed > total_blocks
        chosen_blocks = rng.choice(total_blocks, size=n_blocks_needed,
                                   replace=replace)

        ibd_mean, ibd_var = resample_ibd_with_jackknife_variance(
            pooled_ibd, n_blocks_total=total_blocks,
            pop_ids=pop_ids, pair_info=pair_info, bins=bins,
            target_cm=cm_val, block_size_cm=BLOCK_SIZE_CM,
            chosen_blocks=chosen_blocks,
        )

        leaves = dem_sim.initial_leaves
        n_haploid = np.array([SAMPLES_PER_POP[p] * 2 for p in leaves])
        w_hat, w_se = resample_snp_covariance(
            pooled_snp, chosen_blocks,
            n_haploid_per_pop=n_haploid, se_block_size=50,
        )

        for m_label, m_stan in zip(model_labels, model_stans):
            rep_elbos = [np.nan] * n_topos

            for ti, dem_infer in enumerate(topology_dems):
                n_events_t = len(dem_infer.ordered_events)
                n_admix_t  = dem_infer.n_admix

                init = {"times": [100.0] * n_events_t}
                if n_admix_t > 0:
                    init["admixture_fractions"] = [0.5] * n_admix_t
                if m_label == "IBD-only":
                    init["effective_N"] = float(NE_UNIFORM)
                elif m_label == "Mixed":
                    init["effective_N"] = float(NE_UNIFORM)
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

                    all_vars = fit.stan_variables()
                    cum_t   = all_vars['cumulative_times']
                    admix_f = all_vars.get('admixture_fractions', None)
                    pdict = {}
                    for child, ct_idx, af_idx in admix_mappings[ti]:
                        pdict[f'{child}_time'] = cum_t[:, ct_idx].mean()
                        if admix_f is not None:
                            if admix_f.ndim > 1:
                                pdict[f'{child}_frac'] = admix_f[:, af_idx].mean()
                            else:
                                pdict[f'{child}_frac'] = float(admix_f.mean())
                    if bP_close_idx[ti] is not None:
                        pdict['bP_time'] = cum_t[:, bP_close_idx[ti]].mean()
                    if ab_close_idx[ti] is not None:
                        pdict['ab_time'] = cum_t[:, ab_close_idx[ti]].mean()
                    admix_params[m_label][ti][cm_val].append(pdict)
                except Exception as e:
                    if rep < 3:
                        print(f"    [WARN] {m_label} {topology_labels[ti]} "
                              f"rep {rep+1}: {e}")

            if all(np.isfinite(rep_elbos)):
                elbo_results[m_label][cm_val].append(tuple(rep_elbos))


# ====================================================================
# 9. Plots and summary (verbatim from true-IBD version, only file names change)
# ====================================================================
fig, axes = plt.subplots(2, 3, figsize=(28, 16))
x = np.arange(len(cm_values))
comparisons = [
    ("T_true - T_alt1 (detects RECENT admix?)", 0, 1),
    ("T_true - T_alt2 (detects ANCIENT admix?)", 0, 2),
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
        vp = ax.violinplot(diff_data, positions=x, widths=0.6,
                           showmeans=True, showmedians=True, showextrema=False)
        for body in vp['bodies']:
            body.set_facecolor(color); body.set_alpha(0.5)
            body.set_edgecolor("k"); body.set_linewidth(0.5)
        vp['cmeans'].set_color("k"); vp['cmeans'].set_linewidth(1)
        vp['cmedians'].set_color("k"); vp['cmedians'].set_linewidth(1)
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
    f"Hap-IBD Nfixed (uniform Ne={NE_UNIFORM}): "
    f"T_true vs T_alt1 (no recent) vs T_alt2 (no ancient)\n"
    f"recent on a: t={T_LOOP_SPLIT}->{T_APB_MERGE}->{T_LOOP_CLOSE}, "
    f"alpha={ALPHA_LOOP}  |  ancient on c: t={T_OLD}, f={F_OLD}  |  "
    f"bins from {bins[0][0]} cM",
    fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("admix_thru_b_plus_ancient_dELBO_Nfixed_hapibd.png",
            dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_dELBO_Nfixed_hapibd.png")


# Raw ELBO per topology per model
TOPO_COLORS_ALL = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}
TOPO_NAMES_ALL  = {0: "T_true", 1: "T_alt1", 2: "T_alt2"}
fig_e, axes_e = plt.subplots(1, 3, figsize=(28, 8))
bar_width = 0.8 / n_topos
for m_idx, m_label in enumerate(model_labels):
    ax = axes_e[m_idx]
    for ti in range(n_topos):
        color = TOPO_COLORS_ALL[ti]
        data_list, positions = [], []
        for ci, cm_val in enumerate(cm_values):
            elbos = elbo_results[m_label][cm_val]
            if len(elbos) > 0:
                arr = np.array(elbos); data_list.append(arr[:, ti])
            else:
                data_list.append(np.array([np.nan]))
            positions.append(ci + (ti - (n_topos - 1) / 2) * bar_width)
        vp = ax.violinplot(data_list, positions=positions,
                           widths=bar_width * 0.85,
                           showmeans=True, showmedians=False, showextrema=False)
        for body in vp['bodies']:
            body.set_facecolor(color); body.set_alpha(0.55)
            body.set_edgecolor("k"); body.set_linewidth(0.5)
        vp['cmeans'].set_color(color); vp['cmeans'].set_linewidth(1.5)
    ax.set_xticks(range(len(cm_values)))
    ax.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=10)
    ax.set_xlabel("Genome length (cM)", fontsize=11)
    ax.set_ylabel("ELBO", fontsize=11)
    ax.set_title(f"{m_label}: raw ELBO by topology", fontsize=11, fontweight="bold")
    ax.grid(axis="y", alpha=0.15)
    handles = [Patch(facecolor=TOPO_COLORS_ALL[ti], alpha=0.55,
                     label=TOPO_NAMES_ALL[ti]) for ti in range(n_topos)]
    ax.legend(handles=handles, fontsize=9, loc='best')
fig_e.suptitle(f"Hap-IBD Nfixed (Ne={NE_UNIFORM}): "
               f"Raw ELBO per topology per model",
               fontsize=13, fontweight="bold")
fig_e.tight_layout(rect=[0, 0, 1, 0.94])
fig_e.savefig("admix_thru_b_plus_ancient_ELBO_Nfixed_hapibd.png",
              dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_ELBO_Nfixed_hapibd.png")


# Parameter recovery
TOPO_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}
param_rows = [
    ("Recent admix split time (a)", "a_time",  float(T_LOOP_SPLIT),
     [(0, "T_true"), (2, "T_alt2")]),
    ("aP2 + b merge time (bP)",     "bP_time", float(T_APB_MERGE),
     [(0, "T_true"), (2, "T_alt2")]),
    ("bP + aP1 merge time (ab)",    "ab_time", float(T_LOOP_CLOSE),
     [(0, "T_true"), (2, "T_alt2")]),
    ("Recent admix fraction (a)",   "a_frac",  ALPHA_LOOP,
     [(0, "T_true"), (2, "T_alt2")]),
    ("Ancient admix time (c)",      "c_time",  float(T_OLD),
     [(0, "T_true"), (1, "T_alt1")]),
    ("Ancient admix frac (c)",      "c_frac",  F_OLD,
     [(0, "T_true"), (1, "T_alt1")]),
]
fig2, axes2 = plt.subplots(len(param_rows), 3, figsize=(24, 5 * len(param_rows)))
for row_idx, (row_label, param_key, truth, topo_entries) in enumerate(param_rows):
    for m_idx, m_label in enumerate(model_labels):
        ax = axes2[row_idx, m_idx]
        n_g = len(topo_entries); bar_width = 0.7 / n_g
        for gi, (ti, topo_disp) in enumerate(topo_entries):
            color = TOPO_COLORS[ti]
            data_list, positions = [], []
            for ci, cm_val in enumerate(cm_values):
                vals = [p[param_key] for p in admix_params[m_label][ti][cm_val]
                        if param_key in p]
                data_list.append(np.array(vals) if len(vals) > 0
                                 else np.array([np.nan]))
                positions.append(ci + (gi - (n_g - 1) / 2) * bar_width)
            vp = ax.violinplot(data_list, positions=positions,
                               widths=bar_width * 0.85,
                               showmeans=True, showmedians=False,
                               showextrema=False)
            for body in vp['bodies']:
                body.set_facecolor(color); body.set_alpha(0.55)
                body.set_edgecolor("k"); body.set_linewidth(0.5)
            vp['cmeans'].set_color(color); vp['cmeans'].set_linewidth(1.5)
        ax.axhline(truth, color='red', ls='--', lw=1.2,
                   label=f"truth = {truth:g}")
        ax.set_xticks(range(len(cm_values)))
        ax.set_xticklabels([f"{cm}" for cm in cm_values], fontsize=9)
        ax.set_xlabel("Genome length (cM)", fontsize=10)
        ax.set_ylabel(param_key, fontsize=10)
        ax.set_title(f"{m_label}: {row_label}", fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.15)
        if m_idx == 0:
            handles = [Patch(facecolor=TOPO_COLORS[ti], alpha=0.55, label=name)
                       for ti, name in topo_entries]
            handles.append(plt.Line2D([0], [0], color='red', ls='--',
                                       lw=1.2, label="truth"))
            ax.legend(handles=handles, fontsize=8, loc='best')
fig2.suptitle(f"Hap-IBD Nfixed (Ne={NE_UNIFORM}): parameter recovery",
              fontsize=13, fontweight="bold")
fig2.tight_layout(rect=[0, 0, 1, 0.96])
fig2.savefig("admix_thru_b_plus_ancient_params_Nfixed_hapibd.png",
             dpi=200, bbox_inches="tight")
plt.show()
print("Saved: admix_thru_b_plus_ancient_params_Nfixed_hapibd.png")


# Summary table
print("\n" + "=" * 110)
print("HAP-IBD NFIXED SANITY CHECK - SUMMARY")
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
            d1 = arr[:, 0] - arr[:, 1]; d2 = arr[:, 0] - arr[:, 2]
            pct1 = 100.0 * np.mean(d1 > 0); pct2 = 100.0 * np.mean(d2 > 0)
            print(f"{cm_val:6d} | {m_label:>10s} | "
                  f"{np.mean(d1):+17.1f} | {pct1:9.1f}% | "
                  f"{np.mean(d2):+17.1f} | {pct2:9.1f}% | "
                  f"{len(elbos):5d}")
        else:
            print(f"{cm_val:6d} | {m_label:>10s} | "
                  f"{'N/A':>18s} | {'N/A':>10s} | "
                  f"{'N/A':>18s} | {'N/A':>10s} | {0:5d}")