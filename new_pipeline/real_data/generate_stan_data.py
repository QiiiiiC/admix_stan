"""
Generate Stan-model input data from REAL 1000-Genomes data.

This is the real-data analogue of the simulation pipeline in
``new_pipeline/simulation_methods.py`` / ``bootstrap_ibd.py`` / ``ibd_jackknife.py``,
which all operate on simulated ``tskit`` tree sequences.  Here the same two data
products are produced directly from the files in ``real_data/merged_pruned``:

  1. SNP covariance  (w_hat, w_se)   -- the double-centered TreeMix covariance,
     computed by exactly the same estimator the simulation runs use,
     ``inference_methods.resample_snp_covariance`` (pooled "ratio of sums"
     TreeMix normalization + 1/n_haploid finite-sample bias correction +
     sub-block jackknife SE), fed a per-population allele-frequency matrix read
     from the phased VCFs.

  2. IBD sharing fractions (ibd_mean, ibd_var) per length bin -- IBD segments
     called by hap-IBD / FastSMC are binned by genetic length and turned into
     mean IBD fraction per population pair, with a delete-one-haplotype jackknife
     variance, exactly matching ``ibd_jackknife.resample_ibd_with_jackknife_variance``.

Why this does not reuse the dense (n_pairs x n_blocks) packing of the simulation
pipeline: with 2504 individuals AFR alone has ~870k haplotype pairs, so the dense
arrays are infeasible.  For a *full* real genome we never block-resample to a
target length -- we use every block once -- so the jackknife only needs, per
(bin, pop-pair):  the total summed fraction S, and the per-haplotype sum of
pair-fractions.  Both are accumulated sparsely (one pass over the segments),
and the jackknife formulas are evaluated in closed form including the (huge,
implicit) population of zero-IBD pairs.  This is numerically identical to the
simulation code's jackknife, just without materializing the zeros.

Output (``--out`` prefix):
  <prefix>.npz   : w_hat, w_se, ibd_mean (n_bins,P,P), ibd_var (n_bins,P,P)
  <prefix>.json  : metadata (pop order, haploid counts per pop, bins, genome cM)

The saved matrices are ordered by ``pop_order`` (default AFR,AMR,EAS,EUR,SAS).
To assemble a cmdstanpy data dict, define a ``DemographicTopology`` whose
``initial_leaves`` are in this same order and pass the matrices to
``build_snp_stan_data`` / ``build_ibd_stan_data`` / ``build_mixed_stan_data``
(see ``assemble_stan_data`` helper below).

The IBD caller is chosen with ``--ibd-method``:
  hapibd  -- run hap-IBD (needs the JAR on PATH or HAPIBD_CMD set)
  fastsmc -- run FastSMC via the ``asmc`` Python package
  file    -- read precomputed segment files (with --ibd-glob)
  none    -- SNP product only (default)

Run with the project's stan env, e.g.:
  PY=/opt/anaconda3/envs/stan_env/bin/python

  # SNP only:
  $PY generate_stan_data.py --chrom 22 --out real_stan_data_chr22

  # End-to-end with hap-IBD on chr22:
  $PY generate_stan_data.py --chrom 22 --ibd-method hapibd --out real_chr22

  # End-to-end with FastSMC on chr22:
  $PY generate_stan_data.py --chrom 22 --ibd-method fastsmc --out real_chr22
"""

import os
import sys
import gzip
import glob
import json
import time
import shutil
import argparse
import subprocess
from collections import defaultdict

import numpy as np

# Reuse the new-pipeline SNP covariance estimator so the real-data SNP product
# is computed by exactly the same code path as the simulation run files
# (resample_snp_covariance), not the older per-SNP-standardized variant.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_THIS_DIR)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)
from inference_methods import resample_snp_covariance  # noqa: E402


# Default genome-wide length bins (cM), matching the simulation run files
# exactly (1.0 cM lower edge .. 50 cM) so real- and simulated-data parameter
# recovery use the same bins.
DEFAULT_BINS = [
    [1.0, 1.25], [1.25, 1.5], [1.5, 1.75], [1.75, 2.0], [2.0, 2.5],
    [2.5, 3.0], [3.0, 4.0], [4.0, 5.0], [5.0, 6.0], [6.0, 7.5],
    [7.5, 12.0], [12.0, 20.0], [20.0, 300],
]

DEFAULT_POP_ORDER = ["AFR", "AMR", "EAS", "EUR", "SAS"]

PLOIDY = 2

# Smallest reported segment length (cM): default to the lowest bin boundary so
# every bin can be populated by the IBD caller.
DEFAULT_MIN_CM = min(b[0] for b in DEFAULT_BINS)

# ---- FastSMC (ASMC suite) defaults, mirroring the simulation runner ----
FASTSMC_DEMOGRAPHY     = "CEU"                      # built-in code or path to .demo
FASTSMC_DISCRETIZATION = [(30.0, 100), (100.0, 5)]  # (interval, count) pairs
FASTSMC_DISC_EXTRA     = 39
FASTSMC_TIME           = 100                        # IBD time threshold (generations)
FASTSMC_UKBB_FRQ_URL   = ("https://raw.githubusercontent.com/PalamaraLab/"
                          "PrepareDecoding/main/test/regression/input_UKBB.frq")
# Nominal rates for the SMC++/Oxford map (FastSMC uses the cM column for the
# genetic position; the rate columns are only placeholders for real maps).
FASTSMC_MUT_RATE       = 1.65e-8


# ----------------------------------------------------------------------
# IO helpers
# ----------------------------------------------------------------------
def open_any(path):
    """Open a plain or gzip/bgzip file transparently in text mode."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_population_labels(labels_file, pop_order):
    """
    Read ``population_labels_simple.txt`` (``<sample_id> <POP>`` per line).

    Returns
    -------
    sample_to_pop : dict  {sample_id: pop_idx}
        pop_idx is the index into ``pop_order``.
    """
    pop_to_idx = {p: i for i, p in enumerate(pop_order)}
    sample_to_pop = {}
    with open_any(labels_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            sample_id, pop = parts[0], parts[1]
            if pop not in pop_to_idx:
                raise ValueError(
                    f"Population {pop!r} for sample {sample_id!r} not in "
                    f"pop_order {pop_order}."
                )
            sample_to_pop[sample_id] = pop_to_idx[pop]
    print(f"[labels] {len(sample_to_pop)} samples across {len(pop_order)} pops "
          f"({pop_order})")
    return sample_to_pop


def read_vcf_samples(vcf_path):
    """Return the ordered list of sample IDs from a VCF's #CHROM header line."""
    with open_any(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header found in {vcf_path}")


# ----------------------------------------------------------------------
# Haplotype indexing
# ----------------------------------------------------------------------
def build_haplotype_index(sample_names, sample_to_pop, pop_order):
    """
    Assign every (sample, haplotype) a local index within its population.

    Each diploid sample contributes ``PLOIDY`` haplotypes.  Haplotypes are
    numbered 0..n_pop-1 within each population, ordered by VCF sample order
    then haplotype slot.

    Returns
    -------
    hap_local : dict  {(sample_id, hap_slot): (pop_idx, local_idx)}
        hap_slot is 0-based (hap-IBD's "1"/"2" map to 0/1).
    n_haps : np.ndarray (n_pops,)  haploid counts per population.
    """
    n_pops = len(pop_order)
    counters = [0] * n_pops
    hap_local = {}
    for s in sample_names:
        if s not in sample_to_pop:
            # Sample present in VCF but unlabeled -> skip entirely.
            continue
        p = sample_to_pop[s]
        for h in range(PLOIDY):
            hap_local[(s, h)] = (p, counters[p])
            counters[p] += 1
    n_haps = np.array(counters, dtype=int)
    for p, name in enumerate(pop_order):
        print(f"[haplotypes] {name}: {n_haps[p]} haplotypes "
              f"({n_haps[p] // PLOIDY} individuals)")
    return hap_local, n_haps


# ----------------------------------------------------------------------
# Genetic map
# ----------------------------------------------------------------------
def load_genmap(genmap_dir, chrom):
    """
    Load a PLINK ``.map`` (cols: chrom, snp, cM, bp) for one chromosome.

    Returns (bp_positions, cM_positions) sorted by bp, for np.interp.
    """
    pattern = os.path.join(genmap_dir, f"plink.chr{chrom}.GRCh37.map")
    if not os.path.exists(pattern):
        raise FileNotFoundError(f"Genetic map not found: {pattern}")
    bp, cm = [], []
    with open_any(pattern) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            cm.append(float(parts[2]))
            bp.append(float(parts[3]))
    bp = np.asarray(bp)
    cm = np.asarray(cm)
    order = np.argsort(bp)
    return bp[order], cm[order]


def bp_to_cm(bp_pos, map_bp, map_cm):
    """Interpolate bp -> cM (flat extrapolation outside the mapped range)."""
    return np.interp(bp_pos, map_bp, map_cm)


def chromosome_cm_span(map_bp, map_cm):
    """Total genetic length (cM) of a chromosome from its map."""
    return float(map_cm[-1] - map_cm[0])


# ----------------------------------------------------------------------
# 1. SNP covariance from VCFs
# ----------------------------------------------------------------------
def compute_allele_freq_matrix(vcf_paths, sample_names, sample_to_pop, pop_order,
                               max_sites=None):
    """
    Stream phased VCFs and build a per-population ALT allele-frequency matrix.

    Returns
    -------
    F : np.ndarray (n_snps, n_pops)
        Per-SNP ALT allele frequency in each population.
    """
    n_pops = len(pop_order)

    # Column index lists per population (into the genotype fields parts[9:]).
    pop_cols = [[] for _ in range(n_pops)]
    for col, s in enumerate(sample_names):
        if s in sample_to_pop:
            pop_cols[sample_to_pop[s]].append(col)
    pop_cols = [np.asarray(c, dtype=int) for c in pop_cols]
    # 2 alleles per diploid sample.
    pop_n_alleles = np.array([PLOIDY * len(c) for c in pop_cols], dtype=float)

    freqs_rows = []
    n_seen = 0
    for vcf in vcf_paths:
        with open_any(vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                gts = np.asarray(parts[9:], dtype="U3")
                # Phased diploid GT "a|b"; count ALT alleles fast via char compare.
                # ALT dose per sample = (#'1' among the two slots).
                alt = ((gts == "1|0").astype(np.int8)
                       + (gts == "0|1").astype(np.int8)
                       + 2 * (gts == "1|1").astype(np.int8))
                row = np.empty(n_pops)
                for p in range(n_pops):
                    if pop_n_alleles[p] == 0:
                        row[p] = np.nan
                    else:
                        row[p] = alt[pop_cols[p]].sum() / pop_n_alleles[p]
                freqs_rows.append(row)
                n_seen += 1
                if max_sites is not None and n_seen >= max_sites:
                    break
        if max_sites is not None and n_seen >= max_sites:
            break
        print(f"[snp] processed {vcf} (cumulative {n_seen} sites)")

    F = np.asarray(freqs_rows)
    print(f"[snp] allele-frequency matrix: {F.shape[0]} sites x {F.shape[1]} pops")
    return F


def compute_snp_covariance(vcf_paths, sample_names, sample_to_pop, pop_order,
                           n_haps, se_block_size=500, min_maf=0.0,
                           max_sites=None):
    """Build (w_hat, w_se) with the SAME estimator the simulation runs use.

    This mirrors ``inference_methods.resample_snp_covariance`` exactly -- the
    pooled TreeMix "ratio of sums" covariance ``sum(dev_i dev_j) / sum(het)``,
    the 1/n_haploid finite-sample bias correction, and the ``se_block_size``-SNP
    sub-block jackknife SE -- rather than the per-SNP-standardized
    ``calculate_treemix_covariance``.  Per-SNP deviations are formed against the
    cross-population mean (``mu = mean over pops``), exactly as
    ``inference_methods.prepare_snp_blocks`` does for the simulated data.

    For a full real genome we use every SNP once (no block resampling), so the
    (dev, het) pairs are handed over as a single genome-wide block:
    ``resample_snp_covariance`` re-chunks the concatenation into ``se_block_size``
    sub-blocks for the jackknife SE regardless, so this reproduces the run-file
    estimate exactly.
    """
    F = compute_allele_freq_matrix(
        vcf_paths, sample_names, sample_to_pop, pop_order, max_sites=max_sites
    )

    # Pooled (cross-population) mean per SNP, matching prepare_snp_blocks.
    mu = np.nanmean(F, axis=1)
    valid = (
        (mu > min_maf) & (mu < 1.0 - min_maf)
        & ~np.isnan(mu) & ~np.isnan(F).any(axis=1)
    )
    dev = F[valid] - mu[valid, None]          # (n_snps, n_pops) deviations
    het = mu[valid] * (1.0 - mu[valid])       # (n_snps,) heterozygosities
    print(f"[snp] {int(valid.sum())}/{F.shape[0]} SNPs kept after MAF filter "
          f"(min_maf={min_maf})")

    snp_blocks = [(dev, het)]                  # one genome-wide block
    w_hat, w_se = resample_snp_covariance(
        snp_blocks, chosen_blocks=[0],
        n_haploid_per_pop=np.asarray(n_haps, dtype=float),
        se_block_size=se_block_size,
    )
    print(f"[snp] w_hat range [{w_hat.min():.3e}, {w_hat.max():.3e}], "
          f"w_se range [{w_se.min():.3e}, {w_se.max():.3e}]")
    return w_hat, w_se


# ----------------------------------------------------------------------
# 2. IBD sharing fractions from hap-IBD / FastSMC segments
# ----------------------------------------------------------------------
def _bin_index(length_cm, bins):
    """Return bin index with the convention min < len <= max, else -1."""
    for b_i, (lo, hi) in enumerate(bins):
        if lo < length_cm <= hi:
            return b_i
    return -1


def genmap_file_for(genmap_dir, chrom):
    """Path to the PLINK .map for a chromosome (chrom may be '22' or 'chr22')."""
    chrom = str(chrom)
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return os.path.join(genmap_dir, f"plink.chr{chrom}.GRCh37.map")


# ----------------------------------------------------------------------
# IBD detection -- run the chosen caller on the real VCF
# ----------------------------------------------------------------------
def _run_cmd(cmd, quiet=True):
    if not quiet:
        print(f"  [cmd] {cmd}")
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        print("STDOUT:", res.stdout[-2000:])
        print("STDERR:", res.stderr[-2000:])
        raise RuntimeError(f"Command failed: {cmd}")
    return res


def run_hapibd(vcf_path, map_path, out_prefix, min_seed=DEFAULT_MIN_CM,
               min_output=DEFAULT_MIN_CM, threads=1, quiet=True):
    """
    Run hap-IBD on a real (bgzipped) VCF + PLINK genetic map.

    The command template mirrors the simulation runner.  hap-IBD is a JAR, so
    either put ``hap-ibd`` on PATH or set the ``HAPIBD_CMD`` env var, e.g.:
        export HAPIBD_CMD='java -jar /path/hap-ibd.jar gt={vcf_path} \\
            map={map_path} out={output_prefix} min-seed={min_seed} \\
            min-output={min_output} nthreads={threads}'

    Returns the path to the ``.ibd.gz`` segment file (columns:
    sample1 hap1 sample2 hap2 chrom start_bp end_bp length_cM).
    """
    template = os.environ.get("HAPIBD_CMD")
    if template is None:
        if shutil.which("hap-ibd") is None:
            raise RuntimeError(
                "hap-ibd not found on PATH and HAPIBD_CMD is unset. Install "
                "hap-ibd or set HAPIBD_CMD (see run_hapibd docstring)."
            )
        template = ("hap-ibd gt={vcf_path} map={map_path} out={output_prefix} "
                    "min-seed={min_seed} min-output={min_output} "
                    "nthreads={threads}")

    cmd = template.format(
        vcf_path=vcf_path, map_path=map_path, output_prefix=out_prefix,
        min_seed=min_seed, min_output=min_output, threads=threads,
    )
    t0 = time.time()
    _run_cmd(cmd, quiet=quiet)
    for ext in (".ibd.gz", ".ibd"):
        p = out_prefix + ext
        if os.path.exists(p):
            print(f"  [hap-ibd] {time.time() - t0:.1f}s -> {os.path.basename(p)}")
            return p
    raise FileNotFoundError(f"hap-ibd produced no .ibd[.gz] for {out_prefix}")


def write_oxford_inputs(vcf_path, genmap_dir, chrom, sample_to_pop, out_dir,
                        prefix):
    """
    Convert a real phased VCF into the Oxford inputs FastSMC needs:
      <prefix>.hap.gz  (id1 id2 bp allele0 allele1 <0/1 calls...>)
      <prefix>.samples (ID_1 ID_2 missing header + one row per diploid)
      <prefix>.map     (SMC++ 4-col: bp  rate_cM/Mb  cM  mut_rate)

    Only labeled samples (present in sample_to_pop) are written; their column
    order defines the haplotype order, and the .samples names are echoed back in
    FastSMC's output so downstream (sample, hap) -> pop mapping matches.

    Returns (in_root, sample_names).
    """
    os.makedirs(out_dir, exist_ok=True)
    in_root = os.path.join(out_dir, prefix)
    hap_path = in_root + ".hap.gz"
    smp_path = in_root + ".samples"
    map_path = in_root + ".map"

    sample_names_all = read_vcf_samples(vcf_path)
    keep_cols = [i for i, s in enumerate(sample_names_all) if s in sample_to_pop]
    sample_names = [sample_names_all[i] for i in keep_cols]
    keep_cols = np.asarray(keep_cols, dtype=int)

    map_bp, map_cm = load_genmap(genmap_dir, chrom)

    n_sites = 0
    with gzip.open(hap_path, "wt") as fh:
        with open_any(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                # biallelic SNPs only
                if "," in parts[4]:
                    continue
                bp = int(parts[1])
                gts = np.asarray(parts[9:], dtype="U3")[keep_cols]
                # phased "a|b" -> two 0/1 haplotype calls per diploid
                a = np.where(gts.astype("U1") == "1", "1", "0")
                b = np.array([g[2] if len(g) == 3 else "0" for g in gts])
                b = np.where(b == "1", "1", "0")
                calls = np.empty(2 * len(gts), dtype="U1")
                calls[0::2] = a
                calls[1::2] = b
                fh.write(f"{chrom}:{bp} SNP_{n_sites} {bp} 1 2 "
                         + " ".join(calls) + "\n")
                n_sites += 1

    with open(smp_path, "w") as fs:
        fs.write("ID_1 ID_2 missing\n")
        fs.write("0 0 0\n")
        for nm in sample_names:
            fs.write(f"{nm} {nm} 0\n")

    # SMC++ map: physical bp, cM/Mb rate (placeholder), genetic cM, mut rate.
    with open(map_path, "w") as fm:
        # FastSMC matches by physical position and interpolates the cM column.
        # Provide the genetic-map knots directly (strictly increasing in bp).
        prev_bp = None
        for bp, cm in zip(map_bp.astype(int), map_cm):
            if prev_bp is not None and bp <= prev_bp:
                continue
            fm.write(f"{bp}\t1.0\t{cm:.8f}\t{FASTSMC_MUT_RATE:g}\n")
            prev_bp = bp

    print(f"  [fastsmc] wrote Oxford inputs: {n_sites} sites, "
          f"{len(sample_names)} diploids -> {os.path.basename(in_root)}.*")
    return in_root, sample_names


def build_decoding_quantities(work_root, n_hap, frq_path=None,
                              demography=FASTSMC_DEMOGRAPHY,
                              discretization=FASTSMC_DISCRETIZATION,
                              disc_extra=FASTSMC_DISC_EXTRA):
    """
    Build (and cache) FastSMC decoding quantities for ``n_hap`` haploid samples.

    The CSFS sample size MUST equal the number of haploid samples in the .hap
    file or FastSMC crashes at decode time.  Uses the UKBB frequency *file* so
    any sample size works.  Returns the .decodingQuantities.gz path.
    """
    from asmc.preparedecoding import (
        prepareDecoding, Demography, Discretization, Frequencies,
    )
    os.makedirs(work_root, exist_ok=True)
    dq_root = os.path.join(work_root, f"decoding_quantities_n{n_hap}")
    dq_file = f"{dq_root}.decodingQuantities.gz"
    if os.path.exists(dq_file):
        print(f"  [fastsmc] reusing decoding quantities {os.path.basename(dq_file)}")
        return dq_file

    if frq_path is None:
        frq_path = os.path.join(work_root, "input_UKBB.frq")
    if not os.path.exists(frq_path):
        import urllib.request
        print(f"  [fastsmc] downloading UKBB .frq -> {frq_path}")
        urllib.request.urlretrieve(FASTSMC_UKBB_FRQ_URL, frq_path)

    print(f"  [fastsmc] building decoding quantities (n={n_hap} haploid)...")
    dq = prepareDecoding(
        Demography(demography),
        Discretization(discretization, disc_extra),
        Frequencies(frq_path, n_hap),
        samples=n_hap,
    )
    dq.save_decoding_quantities(dq_root)
    return dq_file


def run_fastsmc(in_root, dq_file, out_root, min_cm=DEFAULT_MIN_CM,
                time_threshold=FASTSMC_TIME):
    """Run FastSMC (array mode + GERMLINE hashing); return its IBD output path."""
    from asmc.asmc import DecodingParams, FastSMC

    params = DecodingParams()
    params.decodingQuantFile         = dq_file
    params.inFileRoot                = in_root
    params.mapFile                   = f"{in_root}.map"
    params.outFileRoot               = out_root
    params.decodingModeString        = "array"
    params.usingCSFS                 = True
    params.batchSize                 = 32
    params.recallThreshold           = 3
    params.min_m                     = min_cm
    params.hashing                   = True
    params.FastSMC                   = True
    params.BIN_OUT                   = True
    params.outputIbdSegmentLength    = True
    params.time                      = time_threshold
    params.noConditionalAgeEstimates = True
    assert params.validateParamsFastSMC()

    t0 = time.time()
    FastSMC(params).run()
    cands = sorted(glob.glob(f"{out_root}*"))
    bibd = [p for p in cands if ".bibd" in os.path.basename(p)]
    out = bibd[0] if bibd else (cands[0] if cands else None)
    if out is None:
        raise FileNotFoundError(f"FastSMC produced no output for {out_root}")
    print(f"  [fastsmc] {time.time() - t0:.1f}s -> {os.path.basename(out)}")
    return out


def fastsmc_to_text(bibd_file, chrom, out_txt):
    """
    Convert FastSMC's binary .bibd.gz output into a hap-IBD-style text file
    (sample1 hap1 sample2 hap2 chrom start end length_cM) so the same parser
    handles both methods.
    """
    from asmc.asmc import BinaryDataReader
    reader = BinaryDataReader(bibd_file)
    n = 0
    with open(out_txt, "w") as out:
        while reader.moreLinesInFile():
            ln = reader.getNextLine()
            start = int(getattr(ln, "ibdStart", 0) or 0)
            end = int(getattr(ln, "ibdEnd", 0) or 0)
            out.write(f"{ln.ind1Id}\t{int(ln.ind1Hap)}\t{ln.ind2Id}\t"
                      f"{int(ln.ind2Hap)}\t{chrom}\t{start}\t{end}\t"
                      f"{float(ln.lengthInCentimorgans):.6f}\n")
            n += 1
    print(f"  [fastsmc] wrote {n} segments -> {os.path.basename(out_txt)}")
    return out_txt


def detect_ibd_for_chrom(method, vcf_path, genmap_dir, chrom, sample_to_pop,
                         work_dir, n_hap, min_cm=DEFAULT_MIN_CM,
                         min_seed=2.0, threads=1, frq_path=None):
    """
    Run the chosen IBD caller on one chromosome's VCF and return the path to a
    segment file readable by ``accumulate_ibd_segments``.
    """
    os.makedirs(work_dir, exist_ok=True)
    if method == "hapibd":
        map_path = genmap_file_for(genmap_dir, chrom)
        out_prefix = os.path.join(work_dir, f"hapibd_chr{chrom}")
        return run_hapibd(vcf_path, map_path, out_prefix, min_seed=min_seed,
                          min_output=min_cm, threads=threads)
    if method == "fastsmc":
        in_root, _ = write_oxford_inputs(
            vcf_path, genmap_dir, chrom, sample_to_pop, work_dir,
            f"fastsmc_chr{chrom}",
        )
        dq_file = build_decoding_quantities(work_dir, n_hap, frq_path=frq_path)
        out_root = os.path.join(work_dir, f"fastsmc_out_chr{chrom}")
        bibd = run_fastsmc(in_root, dq_file, out_root, min_cm=min_cm)
        txt = os.path.join(work_dir, f"fastsmc_chr{chrom}.seg")
        return fastsmc_to_text(bibd, chrom, txt)
    raise ValueError(f"Unknown IBD method {method!r}")


def accumulate_ibd_segments(ibd_files, hap_local, n_haps, pop_order, bins,
                            genmap_dir, genome_cm=None, chrom_filter=None,
                            col_id1=0, col_hap1=1, col_id2=2, col_hap2=3,
                            col_chrom=4, col_start=5, col_end=6, col_cm=7,
                            length_from_map=False):
    """
    One pass over hap-IBD / FastSMC segment files, accumulating the sufficient
    statistics for the haplotype jackknife.

    hap-IBD default columns (0-based):
        0 sample1  1 hap1(1/2)  2 sample2  3 hap2(1/2)
        4 chrom    5 start_bp   6 end_bp   7 length_cM

    IBD fractions are normalized by the genetic length of the chromosome(s)
    actually present in the data (NO 50 cM block: that was only an artifact of
    the 50 cM msprime simulations).  By default the normalizer is the summed
    genetic-map span of every chromosome that appears in the segment files, so
    a single-chromosome run is normalized by that one chromosome's length.

    Parameters
    ----------
    genome_cm : float, optional
        Override for the genetic length normalizer (cM).  If None it is derived
        from the genetic maps of the chromosomes seen in the segments.
    chrom_filter : str or None
        If given (e.g. "22"), only segments on this chromosome are used and the
        normalizer is this chromosome's genetic length.
    length_from_map : bool
        If True, ignore col_cm and compute segment length by interpolating
        start/end bp through the chromosome's genetic map.

    Returns
    -------
    S : dict   {(b_i, i, j): float}            total summed IBD fraction
    hap_sum_i : dict {(b_i, i, j): np.ndarray(n_haps[i])}   per-hap sums, side i
    hap_sum_j : dict {(b_i, i, j): np.ndarray(n_haps[j])}   per-hap sums, side j
        (for within-pop pairs i==j both sides write into the same array, stored
         under hap_sum_i; hap_sum_j[(b,i,i)] is None)
    genome_cm : float   genetic length used as the IBD-fraction normalizer
    """
    n_pops = len(pop_order)
    n_bins = len(bins)

    # Genetic maps (cached) for bp->cM and chromosome length.
    map_cache = {}

    def get_map(chrom):
        chrom = str(chrom)
        chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
        if chrom not in map_cache:
            map_cache[chrom] = load_genmap(genmap_dir, chrom)
        return map_cache[chrom]

    chrom_filter_norm = None
    if chrom_filter is not None:
        chrom_filter_norm = str(chrom_filter)
        chrom_filter_norm = (chrom_filter_norm[3:]
                             if chrom_filter_norm.lower().startswith("chr")
                             else chrom_filter_norm)

    # Accumulate raw segment lengths (cM); divide by genome_cm at the end once
    # we know which chromosomes contributed.
    S = defaultdict(float)
    hap_sum_i = {}
    hap_sum_j = {}
    chroms_seen = set()

    def get_buffers(b_i, i, j):
        key = (b_i, i, j)
        if key not in hap_sum_i:
            hap_sum_i[key] = np.zeros(n_haps[i])
            hap_sum_j[key] = None if i == j else np.zeros(n_haps[j])
        return hap_sum_i[key], hap_sum_j[key]

    n_total = n_kept = n_skip_unmapped = n_skip_bin = n_skip_chrom = 0

    for path in ibd_files:
        with open_any(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                n_total += 1

                chrom = parts[col_chrom]
                chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
                if chrom_filter_norm is not None and chrom != chrom_filter_norm:
                    n_skip_chrom += 1
                    continue

                s1, s2 = parts[col_id1], parts[col_id2]
                # hap-IBD haplotype field is 1/2 -> slot 0/1.
                h1 = int(parts[col_hap1]) - 1
                h2 = int(parts[col_hap2]) - 1

                k1 = (s1, h1)
                k2 = (s2, h2)
                if k1 not in hap_local or k2 not in hap_local:
                    n_skip_unmapped += 1
                    continue

                if length_from_map:
                    mb, mc = get_map(chrom)
                    start_bp = float(parts[col_start])
                    end_bp = float(parts[col_end])
                    length_cm = float(
                        bp_to_cm(end_bp, mb, mc) - bp_to_cm(start_bp, mb, mc)
                    )
                else:
                    length_cm = float(parts[col_cm])

                b_i = _bin_index(length_cm, bins)
                if b_i < 0:
                    n_skip_bin += 1
                    continue

                (pi, ai) = hap_local[k1]
                (pj, aj) = hap_local[k2]
                # Canonical pop-pair order i <= j.
                if pi <= pj:
                    i, j, a, b = pi, pj, ai, aj
                else:
                    i, j, a, b = pj, pi, aj, ai

                chroms_seen.add(chrom)
                S[(b_i, i, j)] += length_cm      # raw cM; normalized below
                buf_i, buf_j = get_buffers(b_i, i, j)
                if i == j:
                    # within-pop: both haplotypes contribute to the same array.
                    buf_i[a] += length_cm
                    buf_i[b] += length_cm
                else:
                    buf_i[a] += length_cm
                    buf_j[b] += length_cm

                n_kept += 1
        print(f"[ibd] parsed {path} (cumulative kept {n_kept})")

    # --- Determine the genetic-length normalizer and convert cM -> fraction ---
    if genome_cm is None:
        if chrom_filter_norm is not None:
            mb, mc = get_map(chrom_filter_norm)
            genome_cm = chromosome_cm_span(mb, mc)
        else:
            genome_cm = 0.0
            for c in sorted(chroms_seen):
                mb, mc = get_map(c)
                genome_cm += chromosome_cm_span(mb, mc)
    print(f"[ibd] chromosomes used: {sorted(chroms_seen)} -> "
          f"normalizer {genome_cm:.2f} cM")

    inv = 1.0 / genome_cm
    for key in S:
        S[key] *= inv
    for key in hap_sum_i:
        hap_sum_i[key] *= inv
        if hap_sum_j[key] is not None:
            hap_sum_j[key] *= inv

    print(f"[ibd] segments: {n_total} read, {n_kept} kept, "
          f"{n_skip_unmapped} unmapped-haplotype, {n_skip_bin} out-of-bin, "
          f"{n_skip_chrom} other-chromosome")
    return dict(S), hap_sum_i, hap_sum_j, genome_cm


def _jackknife_within(S, hap_sums, n):
    """Delete-one-haplotype jackknife variance, within-population (sparse)."""
    P = n * (n - 1) / 2.0
    if P <= 0:
        return 0.0
    mu = S / P
    P_rem = P - (n - 1)
    if P_rem <= 0:
        return 0.0

    # Observed haplotypes (nonzero sum) handled explicitly; the rest share a
    # common pseudovalue because their hap_sum is 0.
    nz = hap_sums[hap_sums != 0.0]
    n_silent = n - nz.size

    # theta_h = n*mu - (n-1)*(S - hap_sums[h]) / P_rem
    leave_nz = (S - nz) / P_rem
    theta_nz = n * mu - (n - 1) * leave_nz
    ss = np.sum((theta_nz - mu) ** 2)

    if n_silent > 0:
        leave0 = S / P_rem
        theta0 = n * mu - (n - 1) * leave0
        ss += n_silent * (theta0 - mu) ** 2

    return ss / (n * (n - 1))


def _jackknife_one_side(S, hap_sums, n_side, n_other, P):
    """One-sided contribution of the between-population jackknife (sparse)."""
    if n_side <= 1 or P <= 0:
        return 0.0
    mu = S / P
    P_rem = P - n_other          # each haplotype on this side is in n_other pairs
    if P_rem <= 0:
        return 0.0

    nz = hap_sums[hap_sums != 0.0]
    n_silent = n_side - nz.size

    leave_nz = (S - nz) / P_rem
    theta_nz = n_side * mu - (n_side - 1) * leave_nz
    ss = np.sum((theta_nz - mu) ** 2)

    if n_silent > 0:
        leave0 = S / P_rem
        theta0 = n_side * mu - (n_side - 1) * leave0
        ss += n_silent * (theta0 - mu) ** 2

    return ss / (n_side * (n_side - 1))


def compute_ibd_matrices(S, hap_sum_i, hap_sum_j, n_haps, bins, se_floor=1e-8):
    """
    Turn the accumulated sufficient statistics into ibd_mean / ibd_var dicts.

    Returns
    -------
    ibd_mean : dict {b_i: np.ndarray (n_pops, n_pops)}
    ibd_var  : dict {b_i: np.ndarray (n_pops, n_pops)}
    """
    n_pops = len(n_haps)
    n_bins = len(bins)

    ibd_mean = {}
    ibd_var = {}
    for b_i in range(n_bins):
        mean_mat = np.zeros((n_pops, n_pops))
        var_mat = np.full((n_pops, n_pops), se_floor ** 2)

        for i in range(n_pops):
            for j in range(i, n_pops):
                n_i = int(n_haps[i])
                n_j = int(n_haps[j])
                if i == j:
                    P = n_i * (n_i - 1) / 2.0
                else:
                    P = float(n_i) * float(n_j)
                if P <= 0:
                    continue

                s = S.get((b_i, i, j), 0.0)
                mean_val = s / P
                mean_mat[i, j] = mean_mat[j, i] = mean_val

                key = (b_i, i, j)
                if i == j:
                    if n_i > 2 and key in hap_sum_i:
                        v = _jackknife_within(s, hap_sum_i[key], n_i)
                    else:
                        v = se_floor ** 2
                else:
                    if n_i > 1 and n_j > 1 and key in hap_sum_i:
                        v_i = _jackknife_one_side(s, hap_sum_i[key], n_i, n_j, P)
                        v_j = _jackknife_one_side(s, hap_sum_j[key], n_j, n_i, P)
                        v = v_i + v_j
                    else:
                        v = se_floor ** 2

                var_mat[i, j] = var_mat[j, i] = max(v, se_floor ** 2)

        ibd_mean[b_i] = mean_mat
        ibd_var[b_i] = var_mat

    return ibd_mean, ibd_var


# ----------------------------------------------------------------------
# Saving / assembling
# ----------------------------------------------------------------------
def save_outputs(out_prefix, pop_order, n_haps, bins, genome_cm,
                 w_hat=None, w_se=None, ibd_mean=None, ibd_var=None):
    """Save matrices to <prefix>.npz and metadata to <prefix>.json."""
    n_bins = len(bins)
    n_pops = len(pop_order)

    arrays = {}
    if w_hat is not None:
        arrays["w_hat"] = w_hat
        arrays["w_se"] = w_se
    if ibd_mean is not None:
        arrays["ibd_mean"] = np.stack([ibd_mean[b] for b in range(n_bins)])
        arrays["ibd_var"] = np.stack([ibd_var[b] for b in range(n_bins)])

    npz_path = out_prefix + ".npz"
    np.savez_compressed(npz_path, **arrays)

    meta = {
        "pop_order": list(pop_order),
        "n_haps": [int(x) for x in n_haps],
        "n_samples_per_pop": {pop_order[p]: int(n_haps[p]) for p in range(n_pops)},
        "bins": [list(b) for b in bins],
        "genome_cm": float(genome_cm) if genome_cm is not None else None,
        "has_snp": w_hat is not None,
        "has_ibd": ibd_mean is not None,
        "arrays_in_npz": list(arrays.keys()),
    }
    json_path = out_prefix + ".json"
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"[save] wrote {npz_path} and {json_path}")
    return npz_path, json_path


def assemble_stan_data(dem, npz_path, json_path, model="mixed", T_max=None,
                       se_floor=1e-8):
    """
    Convenience: load saved matrices and build a cmdstanpy data dict.

    ``dem.initial_leaves`` MUST be in the same order as the saved ``pop_order``.

    model : "ibd" | "snp" | "mixed"
    """
    from simulation_methods import (
        build_ibd_stan_data, build_snp_stan_data, build_mixed_stan_data,
    )
    with open(json_path) as f:
        meta = json.load(f)
    data = np.load(npz_path)

    pop_order = meta["pop_order"]
    if list(dem.initial_leaves) != list(pop_order):
        raise ValueError(
            f"dem.initial_leaves {list(dem.initial_leaves)} must match saved "
            f"pop_order {pop_order} (same order)."
        )

    bins = meta["bins"]
    n_samples = meta["n_haps"]
    genome_cm = meta["genome_cm"]

    if model in ("ibd", "mixed", "snp") and "ibd_mean" in data.files:
        ibd_mean = {b: data["ibd_mean"][b] for b in range(len(bins))}
        ibd_var = {b: data["ibd_var"][b] for b in range(len(bins))}

    if model == "snp":
        return build_snp_stan_data(dem, data["w_hat"], data["w_se"])
    if model == "ibd":
        return build_ibd_stan_data(
            dem, ibd_mean, ibd_var, bins, n_samples_per_pop=n_samples,
            T_max=T_max, se_floor=se_floor, cm=genome_cm,
        )
    if model == "mixed":
        return build_mixed_stan_data(
            dem, ibd_mean, ibd_var, bins, data["w_hat"], data["w_se"],
            n_samples_per_pop=n_samples, T_max=T_max, se_floor=se_floor,
            cm=genome_cm,
        )
    raise ValueError(f"Unknown model {model!r}")


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate Stan input data (SNP covariance + IBD fractions) "
                    "from real 1000G data."
    )
    ap.add_argument("--data-dir", default=os.path.join(_THIS_DIR, "merged_pruned"),
                    help="Directory containing vcf/, genmap/, population_labels_simple.txt")
    ap.add_argument("--labels", default=None,
                    help="Population labels file (default: <data-dir>/population_labels_simple.txt)")
    ap.add_argument("--vcf-glob", default=None,
                    help="Glob for per-chromosome VCFs (default: <data-dir>/vcf/*.vcf.gz)")
    ap.add_argument("--genmap-dir", default=None,
                    help="Genetic map dir (default: <data-dir>/genmap)")
    ap.add_argument("--ibd-method", choices=["hapibd", "fastsmc", "file", "none"],
                    default="none",
                    help="How to obtain IBD segments: run 'hapibd' or 'fastsmc' "
                         "on the VCF, read precomputed segment files ('file', "
                         "with --ibd-glob), or 'none' to skip the IBD product.")
    ap.add_argument("--ibd-glob", default=None,
                    help="With --ibd-method file: glob for precomputed hap-IBD / "
                         "FastSMC segment files.")
    ap.add_argument("--ibd-work-dir", default=None,
                    help="Scratch dir for IBD caller inputs/outputs "
                         "(default: <data-dir>/../<method>_work).")
    ap.add_argument("--out", default=os.path.join(_THIS_DIR, "real_stan_data"),
                    help="Output prefix (writes <prefix>.npz and <prefix>.json)")
    ap.add_argument("--pop-order", nargs="+", default=DEFAULT_POP_ORDER)
    ap.add_argument("--block-size-snps", type=int, default=500,
                    help="SNPs per sub-block for the TreeMix covariance "
                         "jackknife SE (resample_snp_covariance se_block_size).")
    ap.add_argument("--min-maf", type=float, default=0.0)
    ap.add_argument("--max-sites", type=int, default=None,
                    help="Cap on number of SNP sites (debugging).")
    ap.add_argument("--bins-uniform", nargs=3, type=float, default=None,
                    metavar=("LO", "HI", "WIDTH"),
                    help="Use uniform-width length bins from LO to HI cM with "
                         "step WIDTH (e.g. '2 20.5 0.5' for HapNe-style bins). "
                         "Overrides the default unequal bins.")
    ap.add_argument("--min-cm", type=float, default=DEFAULT_MIN_CM,
                    help="Minimum reported IBD segment length (cM).")
    ap.add_argument("--hapibd-min-seed", type=float, default=DEFAULT_MIN_CM,
                    help="hap-IBD min-seed length (cM).")
    ap.add_argument("--threads", type=int, default=1,
                    help="Threads for the IBD caller.")
    ap.add_argument("--frq", default=None,
                    help="UKBB .frq file for FastSMC decoding quantities "
                         "(default: <work-dir>/input_UKBB.frq, downloaded if absent).")
    ap.add_argument("--length-from-map", action="store_true",
                    help="Compute IBD segment length from start/end bp via the "
                         "genetic map instead of trusting the segment's cM column.")
    ap.add_argument("--chrom", default=None,
                    help="Restrict to a single chromosome (e.g. 22). Selects the "
                         "matching VCF and normalizes IBD fractions by that "
                         "chromosome's genetic length. Required when running a "
                         "caller on a single chromosome; omit to loop chr1..22.")
    ap.add_argument("--skip-snp", action="store_true")
    args = ap.parse_args()

    data_dir = args.data_dir
    labels = args.labels or os.path.join(data_dir, "population_labels_simple.txt")
    genmap_dir = args.genmap_dir or os.path.join(data_dir, "genmap")
    pop_order = args.pop_order
    if args.bins_uniform:
        lo, hi, w = args.bins_uniform
        edges = np.arange(lo, hi + 1e-9, w)
        bins = [[float(edges[k]), float(edges[k + 1])] for k in range(len(edges) - 1)]
        print(f"[bins] uniform {lo}-{hi} cM, width {w}: {len(bins)} bins")
    else:
        bins = DEFAULT_BINS

    if args.vcf_glob:
        vcf_glob = args.vcf_glob
    elif args.chrom:
        vcf_glob = os.path.join(
            data_dir, "vcf", f"merged_all.chr{args.chrom}_pruned.vcf.gz"
        )
    else:
        vcf_glob = os.path.join(data_dir, "vcf", "*.vcf.gz")

    vcf_paths = sorted(glob.glob(vcf_glob))
    if not vcf_paths:
        raise FileNotFoundError(f"No VCFs matched {vcf_glob}")

    sample_to_pop = load_population_labels(labels, pop_order)
    sample_names = read_vcf_samples(vcf_paths[0])
    hap_local, n_haps = build_haplotype_index(sample_names, sample_to_pop, pop_order)

    w_hat = w_se = None
    if not args.skip_snp:
        w_hat, w_se = compute_snp_covariance(
            vcf_paths, sample_names, sample_to_pop, pop_order, n_haps,
            se_block_size=args.block_size_snps, min_maf=args.min_maf,
            max_sites=args.max_sites,
        )

    ibd_mean = ibd_var = None
    genome_cm = None

    if args.ibd_method == "none":
        print("[ibd] --ibd-method none; skipping IBD product.")
    else:
        if args.ibd_method == "file":
            if not args.ibd_glob:
                raise ValueError("--ibd-method file requires --ibd-glob.")
            ibd_files = sorted(glob.glob(args.ibd_glob))
            if not ibd_files:
                raise FileNotFoundError(f"No IBD files matched {args.ibd_glob}")
        else:
            # Run the chosen caller per chromosome on the VCF(s).
            work_dir = args.ibd_work_dir or os.path.join(
                data_dir, os.pardir, f"{args.ibd_method}_work"
            )
            n_hap_total = int(n_haps.sum())

            # Map chromosome -> VCF path.
            chrom_to_vcf = {}
            for p in vcf_paths:
                base = os.path.basename(p)
                # e.g. merged_all.chr22_pruned.vcf.gz
                tag = base.split("chr", 1)[-1].split("_")[0] if "chr" in base else None
                if tag:
                    chrom_to_vcf[tag] = p
            if args.chrom:
                c = str(args.chrom)
                c = c[3:] if c.lower().startswith("chr") else c
                if c not in chrom_to_vcf and len(vcf_paths) == 1:
                    chrom_to_vcf = {c: vcf_paths[0]}
                chroms = [c]
            else:
                chroms = sorted(chrom_to_vcf, key=lambda x: int(x) if x.isdigit() else 99)

            ibd_files = []
            for c in chroms:
                if c not in chrom_to_vcf:
                    raise FileNotFoundError(f"No VCF found for chromosome {c}")
                print(f"[ibd] detecting IBD ({args.ibd_method}) on chr{c}...")
                seg = detect_ibd_for_chrom(
                    args.ibd_method, chrom_to_vcf[c], genmap_dir, c,
                    sample_to_pop, work_dir, n_hap_total,
                    min_cm=args.min_cm, min_seed=args.hapibd_min_seed,
                    threads=args.threads, frq_path=args.frq,
                )
                ibd_files.append(seg)

        S, hap_sum_i, hap_sum_j, genome_cm = accumulate_ibd_segments(
            ibd_files, hap_local, n_haps, pop_order, bins, genmap_dir,
            chrom_filter=args.chrom, length_from_map=args.length_from_map,
        )
        ibd_mean, ibd_var = compute_ibd_matrices(
            S, hap_sum_i, hap_sum_j, n_haps, bins,
        )

    save_outputs(args.out, pop_order, n_haps, bins, genome_cm,
                 w_hat=w_hat, w_se=w_se, ibd_mean=ibd_mean, ibd_var=ibd_var)


if __name__ == "__main__":
    main()
