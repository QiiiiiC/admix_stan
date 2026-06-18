"""Nvarying tree ((a,b),c),d): IBD vs SNP vs Mixed, with IBD segments from
FastSMC (the ASMC suite) instead of true MRCA spans. Mirrors the Nfixed
FastSMC variant. Only the parameter relative-error box plot (Mixed model)
is produced; the topology diagram is also saved here.

Requires the asmc package; decoding quantities are built once and cached.
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
# 0b. FastSMC helpers (IBD detection) -- verbatim from the Nfixed variant
# ====================================================================
# --- FastSMC knobs ---
FASTSMC_MIN_CM        = 1.0                    # report segments >= this cM (= smallest bin)
FASTSMC_DEMOGRAPHY    = "CEU"
FASTSMC_DISCRETIZATION = [(30.0, 100), (100.0, 5)]
FASTSMC_DISC_EXTRA    = 39
# CSFS sample size MUST equal the total haploid count (mismatch -> IndexError).
FASTSMC_CSFS_SAMPLES  = sum(N_HAPLOID_PER_POP.values())
FASTSMC_UKBB_FRQ_URL  = ("https://raw.githubusercontent.com/PalamaraLab/"
                         "PrepareDecoding/main/test/regression/input_UKBB.frq")


_DQ_FILE_CACHE = {"path": None}


def write_asmc_inputs(ts, output_dir: Path, prefix: str):
    """Write FastSMC inputs: gzipped Oxford haps (.hap.gz), .samples, .map.

    Diploid individual i is the pair of consecutive sample nodes (2i, 2i+1),
    named tsk_i with FID == IID, so the (IID, hap) -> node mapping matches
    build_sample_to_node below.

    Three format details FastSMC requires (and that differ from a naive
    PLINK-style writer -- getting any of them wrong crashes FastSMC):
      * The haplotypes file must be gzipped: ASMC's hashing pre-pass opens
        <root>.hap.gz.  Columns: `id1 id2 bp allele0 allele1 <0/1 calls...>`.
      * The genetic map is the SMC++/Oxford 4-column format
            physical_pos(bp)  recomb_rate(cM/Mb)  genetic_pos(cM)  mutation_rate
        tab-separated, with the PHYSICAL POSITION IN COLUMN 0 and the genetic
        position (cM) in column 2 -- NOT PLINK `chr id cM bp`.  It is an
        interpolation table (matched to markers by position, not row index),
        must be strictly increasing in physical position, and starts at bp=1.
      * Recombination is uniform here, so cM(bp) = bp * CM_PER_UNIT exactly.

    Returns (in_file_root, sample_names).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    in_root  = output_dir / prefix
    hap_path = output_dir / f"{prefix}.hap.gz"
    smp_path = output_dir / f"{prefix}.samples"
    map_path = output_dir / f"{prefix}.map"

    sample_nodes = list(ts.samples())
    n_hap = len(sample_nodes)
    if n_hap % 2 != 0:
        raise ValueError(f"Expected an even number of sample nodes, got {n_hap}")
    n_dip = n_hap // 2
    sample_names = [f"tsk_{i}" for i in range(n_dip)]

    G = ts.genotype_matrix()                 # (n_sites, n_hap), columns in samples() order
    positions = ts.tables.sites.position     # bp

    # uniform recombination rate expressed in cM/Mb; mutation rate per bp/gen
    rate_cm_per_mb = RECOMB_RATE * 1e6 * 100

    kept_bp = []
    with gzip.open(hap_path, "wt") as fh:
        for s_idx in range(G.shape[0]):
            row = G[s_idx]
            # keep biallelic 0/1 sites only
            if row.max() > 1 or row.min() < 0:
                continue
            pos_bp = int(positions[s_idx])
            gts = " ".join("1" if v else "0" for v in row)
            fh.write(f"1:{pos_bp} SNP_{s_idx} {pos_bp} 1 2 {gts}\n")
            kept_bp.append(pos_bp)

    with open(map_path, "w") as fm:
        # anchor the map at bp=1 / cM=0 unless a marker already sits there
        if not kept_bp or kept_bp[0] > 1:
            fm.write(f"1\t{rate_cm_per_mb:.10f}\t0\t{MUT_RATE:g}\n")
        for bp in kept_bp:
            fm.write(f"{bp}\t{rate_cm_per_mb:.10f}\t{bp * CM_PER_UNIT:.8f}\t{MUT_RATE:g}\n")

    with open(smp_path, "w") as fs:
        fs.write("ID_1 ID_2 missing\n")
        fs.write("0 0 0\n")
        for nm in sample_names:
            fs.write(f"{nm} {nm} 0\n")

    return str(in_root), sample_names

def build_sample_to_node(ts, sample_names):
    """Map (IID, hap_index_1_or_2) -> tskit sample node ID.

    Consecutive sample nodes form the two haplotypes of each diploid, in the
    same order written by write_asmc_inputs.
    """
    sample_nodes = list(ts.samples())
    if len(sample_nodes) != 2 * len(sample_names):
        raise ValueError(
            f"Sample-node count mismatch: {len(sample_nodes)} nodes vs "
            f"{len(sample_names)} diploid samples"
        )
    mapping = {}
    for i, name in enumerate(sample_names):
        mapping[(name, 1)] = sample_nodes[2 * i]
        mapping[(name, 2)] = sample_nodes[2 * i + 1]
    return mapping

def ensure_ukbb_frq(work_root: Path) -> str:
    """Download & cache the UKBB allele-frequency spectrum (.frq) once."""
    frq = work_root / "input_UKBB.frq"
    if not frq.exists():
        import urllib.request
        print(f"[FastSMC] Downloading UKBB .frq -> {frq}")
        urllib.request.urlretrieve(FASTSMC_UKBB_FRQ_URL, str(frq))
    return str(frq)

def build_decoding_quantities(work_root: Path):
    """Build FastSMC decoding quantities once and cache them on disk.

    The CSFS is built for exactly FASTSMC_CSFS_SAMPLES haploid samples (= the
    data's haploid count) using the UKBB frequency *file* (so any sample size
    works, not just the built-in {50,100,200,300}).  The quantities depend only
    on the demography / discretization / frequency spectrum / sample size -- not
    on any particular simulation -- so they are built once and reused for every
    simulation.

    Returns the path to the .decodingQuantities.gz file.
    """
    if _DQ_FILE_CACHE["path"] is not None:
        return _DQ_FILE_CACHE["path"]

    from asmc.preparedecoding import (
        prepareDecoding, Demography, Discretization, Frequencies,
    )

    n = FASTSMC_CSFS_SAMPLES
    dq_root = work_root / f"decoding_quantities_n{n}"
    dq_file = f"{dq_root}.decodingQuantities.gz"

    if not os.path.exists(dq_file):
        frq = ensure_ukbb_frq(work_root)
        print(f"[FastSMC] Building decoding quantities (n={n} haploid) -> {dq_file}")
        dq = prepareDecoding(
            Demography(FASTSMC_DEMOGRAPHY),
            Discretization(FASTSMC_DISCRETIZATION, FASTSMC_DISC_EXTRA),
            Frequencies(frq, n),
            samples=n,
        )
        dq.save_decoding_quantities(str(dq_root))
    else:
        print(f"[FastSMC] Reusing cached decoding quantities {dq_file}")

    _DQ_FILE_CACHE["path"] = dq_file
    return dq_file

def find_fastsmc_output(out_root: str):
    """Locate FastSMC's IBD output.

    With BIN_OUT=True, FastSMC writes <out_root>.<job>.<njobs>.FastSMC.bibd.gz
    (a gzipped binary file that BinaryDataReader opens directly).
    """
    candidates = sorted(glob.glob(f"{out_root}*"))
    if not candidates:
        raise FileNotFoundError(
            f"No FastSMC output found for prefix {out_root}. "
            f"Dir contents: {sorted(glob.glob(os.path.dirname(out_root) + '/*'))}"
        )
    bibd = [p for p in candidates if ".bibd" in os.path.basename(p)]
    if bibd:
        return bibd[0]
    ibd = [p for p in candidates if "ibd" in os.path.basename(p).lower()]
    return ibd[0] if ibd else candidates[0]

def run_fastsmc(in_root: str, dq_file: str, out_root: str):
    """Run FastSMC (array mode + GERMLINE hashing) and return its IBD output.

    Uses the default-constructed DecodingParams with every attribute set
    explicitly (the 4-arg "FastSMC=True" constructor overload does NOT actually
    flip params.FastSMC on this build, so we set it ourselves), then the
    required validateParamsFastSMC() pre-run step.
    """
    from asmc.asmc import DecodingParams, FastSMC

    params = DecodingParams()
    params.decodingQuantFile      = dq_file
    params.inFileRoot             = in_root
    params.mapFile                = f"{in_root}.map"
    params.outFileRoot            = out_root
    params.decodingModeString     = "array"
    params.usingCSFS              = True
    params.batchSize              = 32
    params.recallThreshold        = 3
    params.min_m                  = FASTSMC_MIN_CM   # min reported segment length (cM)
    params.hashing                = True             # GERMLINE-style hashing pre-pass
    params.FastSMC                = True
    params.BIN_OUT                = True             # binary .bibd.gz output
    params.outputIbdSegmentLength = True             # populate lengthInCentimorgans
    params.time                   = 100              # IBD time threshold (generations)
    params.noConditionalAgeEstimates = True
    assert params.validateParamsFastSMC()

    fastsmc = FastSMC(params)
    fastsmc.run()
    return find_fastsmc_output(out_root)

def parse_fastsmc_segments(output_file: str, sample_to_node):
    """Read FastSMC IBD output into [(node_u, node_v, length_cm), ...].

    Tries the binary reader first (asmc.asmc.BinaryDataReader); falls back to a
    whitespace-delimited text parse with the same column meaning as Hap-IBD
    (id1 hap1 id2 hap2 chrom start end length_cm).
    """
    segments = []

    # --- binary path ---
    try:
        from asmc.asmc import BinaryDataReader
        reader = BinaryDataReader(output_file)
        while reader.moreLinesInFile():
            line = reader.getNextLine()
            key1 = (str(line.ind1Id), int(line.ind1Hap))
            key2 = (str(line.ind2Id), int(line.ind2Hap))
            if key1 in sample_to_node and key2 in sample_to_node:
                segments.append((sample_to_node[key1], sample_to_node[key2],
                                 float(line.lengthInCentimorgans)))
        return segments
    except Exception as e:
        print(f"[FastSMC] binary read failed ({e}); trying text parse")

    # --- text fallback ---
    opener = gzip.open if str(output_file).endswith(".gz") else open
    with opener(output_file, "rt") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            t = ln.split()
            if len(t) < 8:
                continue
            try:
                length_cm = float(t[7])
            except ValueError:
                continue
            key1 = (t[0], int(t[1]))
            key2 = (t[2], int(t[3]))
            if key1 in sample_to_node and key2 in sample_to_node:
                segments.append((sample_to_node[key1], sample_to_node[key2],
                                 length_cm))
    return segments

def fastsmc_segments_to_packed(ts, segments, bins, block_size_cm, cm_per_unit):
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

    pair_info = {}
    pair_key_to_index = {}
    packed_shape = {}

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

    packed = {key: np.zeros((packed_shape[key], n_blocks))
              for key in packed_shape}

    node_to_pop_idx = {sample_nodes[k]: np.searchsorted(pop_ids, node_to_pop[k])
                       for k in range(len(sample_nodes))}

    n_segs_kept = 0
    n_segs_skipped = 0
    for u, v, length_cm in segments:
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

work_root = Path(_PARENT) / "fastsmc_work"
work_root.mkdir(exist_ok=True)
# Per-script tag: isolates each script's per-simulation dirs under the
# shared work root, so two different run scripts can run concurrently
# without clobbering each other's sim_* inputs/outputs.
_RUN_TAG = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]
dq_file = build_decoding_quantities(work_root)
print(f"FastSMC decoding quantities: {dq_file}")

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
    in_root, sample_names = write_asmc_inputs(mts, sim_dir, "sim")
    out_root = str(sim_dir / "fastsmc_out")
    fastsmc_out = run_fastsmc(in_root, dq_file, out_root)
    sample_to_node = build_sample_to_node(mts, sample_names)
    segments = parse_fastsmc_segments(fastsmc_out, sample_to_node)
    packed, n_blocks, pop_samples, pids, pi, n_kept, n_skip = fastsmc_segments_to_packed(
        mts, segments, bins, block_size_cm=BLOCK_SIZE_CM, cm_per_unit=CM_PER_UNIT,
    )
    print(f"  [FastSMC] segs kept={n_kept} skipped={n_skip}")
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
        outpath=os.path.join(_HERE, f"relative_error_fastsmc_{_m_tag}.png"),
        title=f"Parameter relative error ({_m_label}, FastSMC)",
    )
