"""IBD segments from a tree sequence, binned by length into [bin, i, j] arrays.

Two sources, same output shape so everything downstream is source-agnostic:
  true_ibd  strict MRCA-node spans read off the genealogy (C++ scanner)
  hapibd    segments called by hap-IBD from the simulated genotypes

Output: count[b, i, j] (segments) and length[b, i, j] (summed cM), symmetric in
(i, j), for populations in `pops` order, over bins edges[b] < L <= edges[b+1].
Pairs are HAPLOTYPE pairs: within a population n(n-1)/2, between n_i * n_j.
"""
from __future__ import annotations

import ctypes
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
JAVA = os.environ.get("JAVA", "/opt/miniconda3/envs/genetics_env/bin/java")
HAPIBD_JAR = os.environ.get(
    "HAPIBD_JAR", "/opt/miniconda3/envs/genetics_env/share/hap-ibd-1.0.rev20May22.818-0/hap-ibd.jar")


def bin_edges(lo, hi, step):
    """Uniform edges lo, lo+step, ..., hi (cM)."""
    n = (hi - lo) / step
    if not np.isclose(n, round(n)):
        raise ValueError("bin range must be a whole number of steps")
    return np.linspace(lo, hi, round(n) + 1)


def pair_counts(n_haploid):
    """Haplotype-pair exposure per population pair: n(n-1)/2 on the diagonal, n_i n_j off it."""
    n = np.asarray(n_haploid, float)
    out = np.outer(n, n)
    out[np.diag_indices(len(n))] = n * (n - 1) / 2
    return out


def symmetrize(values):
    """Copy the upper triangle of [..., i, j] onto the lower one."""
    i, j = np.triu_indices(values.shape[-1], 1)
    values[..., j, i] = values[..., i, j]
    return values


def sample_pairs(ts, pops):
    """All sample pairs, each pair's upper-triangle cell i*L + j, and sample -> population index."""
    pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
    lookup = {int(u): i for i, name in enumerate(pops) for u in ts.samples(population=pop_id[name])}
    L = len(pops)
    uv = np.asarray(list(itertools.combinations(ts.samples(), 2)), dtype=np.int32)
    cells = np.asarray([min(lookup[int(u)], lookup[int(v)]) * L + max(lookup[int(u)], lookup[int(v)])
                        for u, v in uv], dtype=np.int32)
    return uv, cells, lookup


def mrca_library():
    """Compile (once per source hash) and load the C++ MRCA-span scanner."""
    folder = HERE / ".build"
    folder.mkdir(exist_ok=True)
    src = HERE / "mrca_scan.cpp"
    lib = folder / ("mrca_" + hashlib.sha256(src.read_bytes()).hexdigest()[:12] + ".so")
    if not lib.exists():
        compiler = shutil.which(os.environ.get("CXX", "c++"))
        if not compiler:
            raise RuntimeError("C++ compiler required for fast strict-MRCA extraction (CXX or c++).")
        tmp = lib.with_suffix(f".{os.getpid()}.tmp.so")
        subprocess.run([compiler, "-O3", "-std=c++11", "-shared", "-fPIC", str(src), "-o", str(tmp)], check=True)
        tmp.replace(lib)
    fn = ctypes.CDLL(str(lib)).scan
    ip = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
    dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    lp = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
    fn.argtypes = [ctypes.c_int, ip, ip, ip, ip, dp, ctypes.c_double, ctypes.c_double,
                   ctypes.c_int, dp, ip, dp, lp, dp, ctypes.c_int, ctypes.c_int]
    fn.restype = None
    return fn


def true_ibd(ts, pops, edges, recombination_rate):
    """Segments = maximal stretches over which a pair's MRCA node is unchanged."""
    uv, cells, _ = sample_pairs(ts, pops)
    L = len(pops)
    u, v = (np.ascontiguousarray(uv[:, k]) for k in (0, 1))
    edges = np.ascontiguousarray(edges, dtype=float)
    previous = np.full(len(uv), -2, np.int32)
    start = np.zeros(len(uv))
    count = np.zeros((len(edges) - 1, L, L), np.int64)
    length = np.zeros_like(count, dtype=float)
    fn = mrca_library()
    times = np.ascontiguousarray(ts.nodes_time)
    cm_per_bp = 100 * recombination_rate
    for tree in ts.trees():
        fn(len(u), u, v, cells, tree.parent_array, times, tree.interval.left,
           cm_per_bp, len(edges) - 1, edges, previous, start, count, length, L * L, 0)
    fn(len(u), u, v, cells, tree.parent_array, times, ts.sequence_length,
       cm_per_bp, len(edges) - 1, edges, previous, start, count, length, L * L, 1)
    return symmetrize(count), symmetrize(length)


def hapibd_command():
    """argv prefix for hap-IBD: $HAPIBD_COMMAND, else the genetics_env Java + jar."""
    import shlex
    if os.environ.get("HAPIBD_COMMAND"):
        return shlex.split(os.environ["HAPIBD_COMMAND"])
    return [JAVA, "-jar", HAPIBD_JAR]


def hapibd_version(command):
    p = subprocess.run(command, capture_output=True, text=True, timeout=30)
    text = p.stdout + p.stderr
    if "hap-ibd" not in text.lower() or "Unable to locate a Java Runtime" in text:
        raise RuntimeError(f"Could not launch hap-IBD: {text[-2000:]}")
    return text


def hapibd(ts, pops, edges, recombination_rate, folder, command,
           min_seed=1.0, min_output=2.0, threads=1):
    """Call IBD with hap-IBD on the simulated genotypes (uniform genetic map).

    command is an argv prefix, e.g. ['java', '-jar', '/path/hap-ibd.jar'].
    min_output must not exceed the first bin edge, or the lowest bin is truncated.
    """
    if min_output > edges[0]:
        raise ValueError("hap-IBD min-output above the lowest bin edge")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    vcf, gmap, out = folder / "input.vcf", folder / "input.map", folder / "hap"
    individuals = [ind for ind in ts.individuals() if len(ind.nodes) == 2 and
                   all(ts.node(int(u)).is_sample() for u in ind.nodes)]
    names = [f"sample_{ind.id}" for ind in individuals]
    mapping = {(name, k + 1): int(u) for name, ind in zip(names, individuals)
               for k, u in enumerate(ind.nodes)}
    if len(mapping) != ts.num_samples:
        raise ValueError("Every sample must belong to a diploid individual.")
    with vcf.open("w") as handle:
        ts.write_vcf(handle, contig_id="1", individuals=[x.id for x in individuals],
                     individual_names=names, position_transform=lambda x: np.asarray(x) + 1)
    scale = 100 * recombination_rate
    gmap.write_text(f"1 . 0 1\n1 . {ts.sequence_length*scale:.12g} {int(ts.sequence_length)+1}\n")
    args = list(command) + [f"gt={vcf}", f"map={gmap}", f"out={out}", f"min-seed={min_seed}",
                            f"min-output={min_output}", f"nthreads={threads}"]
    with (folder / "command.log").open("w") as handle:
        handle.write(json.dumps(args) + "\n")
        handle.flush()
        subprocess.run(args, stdout=handle, stderr=subprocess.STDOUT, check=True)
    _, _, lookup = sample_pairs(ts, pops)
    L = len(pops)
    count = np.zeros((len(edges) - 1, L, L), np.int64)
    length = np.zeros_like(count, dtype=float)
    # .hbd.gz holds the homologous pair inside each diploid; the diagonal
    # exposure n(n-1)/2 includes those pairs, so BOTH files are required.
    for suffix in (".ibd.gz", ".hbd.gz"):
        with gzip.open(str(out) + suffix, "rt") as handle:
            for line in handle:
                x = line.split()
                if len(x) != 8:
                    raise ValueError(f"Malformed hap-IBD output: {line!r}")
                u = mapping[(x[0], int(x[1]))]
                v = mapping[(x[2], int(x[3]))]
                if u == v:
                    raise ValueError("Unexpected same-haplotype pair")
                i, j = sorted((lookup[u], lookup[v]))
                size = float(x[7])
                b = np.searchsorted(edges, size, side="left") - 1
                if 0 <= b < len(count):
                    count[b, i, j] += 1
                    length[b, i, j] += size
    return symmetrize(count), symmetrize(length)
