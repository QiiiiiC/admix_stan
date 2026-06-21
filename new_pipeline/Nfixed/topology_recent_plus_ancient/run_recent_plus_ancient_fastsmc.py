"""
FastSMC version of the Nfixed recent + ancient admixture experiment.

Mirrors run_recent_plus_ancient_hapibd.py exactly, except IBD segments
come from FastSMC (the ASMC suite, via the `asmc` Python package) instead of
Hap-IBD.

Pipeline (per simulation, identical downstream to the Hap-IBD version):
  1. msprime simulate 50 cM at human-realistic params -> tree sequence
  2. Write Oxford phased input (.hap/.samples) + a PLINK genetic map
  3. Run FastSMC (asmc.asmc.FastSMC) using precomputed decoding quantities
  4. Parse the binary IBD output, map (IID, hap) -> tskit haplotype node, and
     pack into packed[(b_i, i, j)] in the same deterministic pair order as
     calculate_ibd_blocks_mrca / the Hap-IBD version.
  5. SNP blocks computed exactly as before.

After all sims, pool with pool_multiple_simulations and resample with
resample_ibd_with_jackknife_variance --- both unchanged.

Only three figures are produced (shared with the true-IBD and Hap-IBD scripts):
  1. dELBO comparison box plot (all 3 models, with means)        [admix_plots]
  2. parameter comparison, IBD-only vs Mixed                     [admix_plots]
  3. parameter relative error, Mixed model only                  [relative_error]

Requirements:
  - asmc Python package (`pip install asmc-asmc asmc-preparedecoding`)
  - Decoding quantities are built once via asmc.preparedecoding.prepareDecoding
    and cached on disk; tweak the FASTSMC_* constants below to match the
    demography / discretization / frequencies you want.

NOTE on FastSMC I/O: FastSMC's exact input-file expectations and the name of
its output file can vary between asmc builds. The detection is isolated in
write_asmc_inputs() / build_decoding_quantities() / run_fastsmc() /
parse_fastsmc_segments(); adjust those if your asmc build differs. Everything
downstream is shared, verified code.
"""

import sys, os, time, gzip, shutil, glob, re
from pathlib import Path
from collections import defaultdict

_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _PARENT)
_MODELS = os.path.join(_PARENT, "models")
# Write every figure (results + topology diagram) next to this script.
_HERE = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import matplotlib.pyplot as plt

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
from relative_error import plot_relative_error_boxplot
from admix_plots import plot_delbo_boxplot, plot_param_comparison
from cmdstanpy import CmdStanModel


# ====================================================================
# 0. Configuration  (human-realistic units, identical to the Hap-IBD run)
# ====================================================================
N_REPLICATES   = 100
BLOCK_SIZE_CM  = 50.0
CM_PER_UNIT    = 1e-6       # 1 cM per Mb (standard human)
RECOMB_RATE    = 1e-8       # crossovers/bp/generation
MUT_RATE       = 1.25e-8    # mutations/bp/generation
SAMPLES_PER_POP = {'a': 20, 'b': 20, 'c': 20, 'd': 20, 'e': 20}
# Haploid sample count per population (diploid samples x ploidy 2)
N_HAPLOID_PER_POP = {_p: 2 * _c for _p, _c in SAMPLES_PER_POP.items()}

N_SIMS = 50
SIM_CM_EACH = 50
SIM_SEQ_LEN = SIM_CM_EACH / CM_PER_UNIT

cm_values = [50, 100, 150, 200, 300, 500, 750, 1000]

# Compare bins from 1.0 cM (matches the Hap-IBD run exactly)
bins = [
    [1.0, 1.25],[1.25,1.5],
    [1.5, 1.75],[1.75,2.0], [2.0, 2.5],[2.5,3.0],
    [3.0,4.0],[4.0,5.0],[5.0, 6.0],[6.0,7.5],
    [7.5, 12.0],[12.0,20.0], [20.0, BLOCK_SIZE_CM]
]

# Single uniform Ne everywhere (haploid; demography API takes diploid = haploid//2)
NE_UNIFORM = 20000

# --- Recent admixture on b ---
T_RECENT     = 30      # b admixture time
ALPHA_RECENT = 0.70    # b: fraction from the a-side parent (bP1)

# --- Ancient admixture on d ---
T_OLD = 500            # d admixture time
F_OLD = 0.60           # d: fraction from the c-side parent (dP1)

# --- Merge event times ---
T_AB    = 50
T_CB    = 100
T_CBD   = 700
T_LEFT  = 1000
T_RIGHT = 1200
T_ROOT  = 1500

SNP_CUTOFF_TIME = T_ROOT
SNP_MIN_MAF = 0.05

# --- FastSMC knobs ---
FASTSMC_MIN_CM        = 1.0                    # report segments >= this cM
FASTSMC_DEMOGRAPHY    = "CEU"                  # built-in code, or path to a .demo file
FASTSMC_DISCRETIZATION = [(30.0, 100), (100.0, 5)]  # ASMC (interval, count) pairs
FASTSMC_DISC_EXTRA    = 39                     # extra quantile-based discretization points
# The CSFS sample size MUST equal the total number of haploid samples in the data.
# ASMC indexes per-site arrays by this count, so a mismatch crashes the decoder
# with "IndexError: vector" at run time (this was the original bug).
FASTSMC_CSFS_SAMPLES  = sum(N_HAPLOID_PER_POP.values())   # 5 pops x 15 dip x 2 = 150
# UKBB allele-frequency spectrum file.  The built-in "UKBB" string only ships
# precomputed spectra for sample sizes {50, 100, 200, 300}; passing the .frq file
# instead lets prepareDecoding build the CSFS for ANY sample size while keeping
# UKBB's (well-calibrated, array-like) frequency spectrum -- which is what gives
# good IBD detection.  A data-derived .frq (sequence SFS) detects almost nothing.
FASTSMC_UKBB_FRQ_URL  = ("https://raw.githubusercontent.com/PalamaraLab/"
                         "PrepareDecoding/main/test/regression/input_UKBB.frq")


# ====================================================================
# 1. FastSMC helpers
# ====================================================================
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
            # +1: positions are 1-based; a site at bp 0 trips FastSMC's
            # "increasing physical position" check (0 vs the initial 0).
            pos_bp = int(positions[s_idx]) + 1
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


_DQ_FILE_CACHE = {"path": None}


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


# ====================================================================
# 2. Pack FastSMC segments into the calculate_ibd_blocks_mrca format
#    (identical to the Hap-IBD packer -- detection-method agnostic)
# ====================================================================
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
# 3. Stan helpers (unchanged)
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
# 4. Topology builders (verbatim from the Hap-IBD script)
# ====================================================================
def _set_uniform_ne(dem, ne_haploid):
    for n in dem.nodes:
        dem.set_node_ne(n, ne_haploid // 2)


def build_t_true():
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    d.add_admixture_event('b', 'bP1', 'bP2')
    d.add_merge_event('a', 'bP1', 'ab')
    d.add_merge_event('c', 'bP2', 'cb')
    d.add_admixture_event('d', 'dP1', 'dP2')
    d.add_merge_event('cb', 'dP1', 'cbd')
    d.add_merge_event('ab', 'cbd', 'left')
    d.add_merge_event('e', 'dP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.set_admixture_parameters('b', time=T_RECENT,
                               fraction_parent_1=ALPHA_RECENT, parent_1_name='bP1')
    d.set_admixture_parameters('d', time=T_OLD,
                               fraction_parent_1=F_OLD, parent_1_name='dP1')
    d.set_merge_time('ab', T_AB)
    d.set_merge_time('cb', T_CB)
    d.set_merge_time('cbd', T_CBD)
    d.set_merge_time('left', T_LEFT)
    d.set_merge_time('right', T_RIGHT)
    d.set_merge_time('root', T_ROOT)
    d.finalize_root()
    return d


def build_t_alt1():
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    d.add_merge_event('a', 'b', 'ab')
    d.add_admixture_event('d', 'dP1', 'dP2')
    d.add_merge_event('c', 'dP1', 'cdP1')
    d.add_merge_event('ab', 'cdP1', 'left')
    d.add_merge_event('e', 'dP2', 'right')
    d.add_merge_event('left', 'right', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.finalize_root()
    return d


def build_t_alt2():
    d = DemographicTopology(['a', 'b', 'c', 'd', 'e'])
    d.add_admixture_event('b', 'bP1', 'bP2')
    d.add_merge_event('a', 'bP1', 'ab')
    d.add_merge_event('c', 'bP2', 'cb')
    d.add_merge_event('cb', 'd', 'cbd')
    d.add_merge_event('ab', 'cbd', 'left')
    d.add_merge_event('left', 'e', 'root')
    _set_uniform_ne(d, NE_UNIFORM)
    d.finalize_root()
    return d


# ====================================================================
# 5. Main pipeline
# ====================================================================
print("=" * 60)
print("FastSMC estimation pipeline -- Nfixed recent + ancient admixture")
print(f"  Sim length: {SIM_CM_EACH} cM  ({SIM_SEQ_LEN} bp)")
print(f"  N sims    : {N_SIMS}")
print(f"  CM per bp : {CM_PER_UNIT}  (recomb {RECOMB_RATE}/bp/gen)")
print("=" * 60)

dem_true = build_t_true()

# --- Save a diagram of the (true) data-generating topology in this folder ---
_topo_fig, _topo_ax = plt.subplots(figsize=(10, 6))
dem_true.plot_demography(scale=True, ax=_topo_ax)
_topo_ax.set_title('recent_plus_ancient (true topology)')
_topo_fig.savefig(os.path.join(_HERE, 'topology_recent_plus_ancient.png'),
                  dpi=200, bbox_inches='tight')
plt.close(_topo_fig)
print(f"[topology] saved {os.path.join(_HERE, 'topology_recent_plus_ancient.png')}")
dem_alt1 = build_t_alt1()
dem_alt2 = build_t_alt2()
topology_dems = [dem_true, dem_alt1, dem_alt2]
topology_labels = [
    "T_true (recent b + ancient d)",
    "T_alt1 (ancient only)",
    "T_alt2 (recent only)",
]
for lbl, d in zip(topology_labels, topology_dems):
    print(f"  {lbl}: n_events={len(d.ordered_events)}, "
          f"n_admix={d.n_admix}, n_nodes={len(d.nodes)}")
n_topos = len(topology_dems)
dem_sim = dem_true


# ---- Per-simulation: msprime -> ASMC inputs -> FastSMC -> packed ----
packed_list     = []
n_blocks_list   = []
snp_blocks_list = []
pair_info       = None
pop_ids         = None

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

    sim_dir = work_root / f"{_RUN_TAG}_sim_{sim_i}"
    in_root, sample_names = write_asmc_inputs(mts, sim_dir, "sim")

    out_root = str(sim_dir / "fastsmc_out")
    fastsmc_out = run_fastsmc(in_root, dq_file, out_root)
    t2 = time.time()
    print(f"  [FastSMC] {t2 - t1:.2f}s -> {os.path.basename(fastsmc_out)}")

    sample_to_node = build_sample_to_node(mts, sample_names)
    segments = parse_fastsmc_segments(fastsmc_out, sample_to_node)

    packed, n_blocks, pop_samples, pids, pi, n_kept, n_skip = \
        fastsmc_segments_to_packed(
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

    # SNP blocks: exactly as in the true-IBD / Hap-IBD versions
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
for ti, (lbl, mp) in enumerate(zip(topology_labels, admix_mappings)):
    print(f"  {lbl} admix mapping: {mp}")


# ====================================================================
# 8. Replicates: resample blocks -> fit Stan for each (cm, rep, topo, model)
# ====================================================================
elbo_results = {label: {cm: [] for cm in cm_values} for label in model_labels}
admix_params = {
    label: {ti: {cm: [] for cm in cm_values} for ti in range(n_topos)}
    for label in model_labels
}
# Per-replicate stan_variables() for the relative-error box plot:
# Mixed model, true topology (ti == 0), largest genome length only.
rel_err_stanvars = []
REL_ERR_CM = cm_values[-1]

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
                # All three Nfixed models carry a single `effective_N`; it must be
                # initialised near the truth or Pathfinder fails from its random
                # start (this was the SNP-only "no result" bug -- SNP-only was the
                # only model not getting an effective_N init).
                init["effective_N"] = float(NE_UNIFORM)
                if m_label == "Mixed":
                    init["kappa_snp"] = 1.0

                try:
                    if m_label == "IBD-only":
                        sd = build_ibd_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, n_samples_per_pop=N_HAPLOID_PER_POP,
                            T_max=100000, cm=cm_val,
                        )
                    elif m_label == "SNP-only":
                        sd = build_snp_stan_data(
                            dem_infer, w_hat, w_se,
                        )
                    else:
                        sd = build_mixed_stan_data(
                            dem_infer, ibd_mean, ibd_var, bins, w_hat, w_se, n_samples_per_pop=N_HAPLOID_PER_POP,
                            T_max=100000, cm=cm_val,
                        )

                    fit = m_stan.pathfinder(
                        data=sd, inits=init, show_console=False,
                        psis_resample=True,
                    )
                    rep_elbos[ti] = extract_elbo(fit)

                    all_vars = fit.stan_variables()

                    # Collect for the relative-error box plot.
                    if (m_label == "Mixed" and ti == 0
                            and cm_val == REL_ERR_CM):
                        rel_err_stanvars.append(all_vars)

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
                    admix_params[m_label][ti][cm_val].append(pdict)
                except Exception as e:
                    if rep < 3:
                        print(f"    [WARN] {m_label} {topology_labels[ti]} "
                              f"rep {rep+1}: {e}")

            if all(np.isfinite(rep_elbos)):
                elbo_results[m_label][cm_val].append(tuple(rep_elbos))


# ====================================================================
# 9. The three figures
# ====================================================================
comparisons = [
    ("T_true - T_alt1 (detects RECENT admix?)", 0, 1),
    ("T_true - T_alt2 (detects ANCIENT admix?)", 0, 2),
]
plot_delbo_boxplot(
    elbo_results, cm_values, model_labels, model_colors, comparisons,
    outpath=os.path.join(_HERE, "recent_plus_ancient_dELBO_Nfixed_fastsmc.png"),
    suptitle=(f"FastSMC Nfixed (uniform Ne={NE_UNIFORM}): dELBO comparison  |  "
              f"recent b admix t={T_RECENT}, alpha={ALPHA_RECENT}  |  "
              f"ancient d admix t={T_OLD}, f={F_OLD}  |  "
              f"bins from {bins[0][0]} cM"),
)

TOPO_COLORS = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}
param_rows = [
    ("Recent admix time (b)",      "b_time", float(T_RECENT),
     [(0, "T_true"), (2, "T_alt2")]),
    ("Recent admix fraction (b)",  "b_frac", ALPHA_RECENT,
     [(0, "T_true"), (2, "T_alt2")]),
    ("Ancient admix time (d)",     "d_time", float(T_OLD),
     [(0, "T_true"), (1, "T_alt1")]),
    ("Ancient admix fraction (d)", "d_frac", F_OLD,
     [(0, "T_true"), (1, "T_alt1")]),
]
plot_param_comparison(
    admix_params, cm_values, param_rows, ["IBD-only", "Mixed"], TOPO_COLORS,
    outpath=os.path.join(_HERE, "recent_plus_ancient_params_Nfixed_fastsmc.png"),
    suptitle=f"FastSMC Nfixed (Ne={NE_UNIFORM}): parameter recovery (IBD-only vs Mixed)",
)

plot_relative_error_boxplot(
    dem_true, rel_err_stanvars, varying=False,
    outpath=os.path.join(_HERE, "recent_plus_ancient_relative_error_Nfixed_fastsmc.png"),
    title=f"Parameter relative error (Mixed, true topology, {REL_ERR_CM} cM)",
)


# ====================================================================
# 10. Summary table
# ====================================================================
print("\n" + "=" * 110)
print("FASTSMC NFIXED SANITY CHECK - SUMMARY")
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
