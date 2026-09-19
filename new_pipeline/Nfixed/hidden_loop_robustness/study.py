"""Shared design, strict MRCA counts, hap-IBD, and matched SNP summaries."""
from __future__ import annotations

import contextlib
import ctypes
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".build" / "matplotlib"))

import numpy as np

HERE = Path(__file__).resolve().parent
PIPELINE = HERE.parents[1]
REAL_DATA = PIPELINE / "real_data"
sys.path[:0] = [str(PIPELINE), str(REAL_DATA / "3pop")]
from demography import DemographicTopology
from enumerate_3pop import enumerate_all, make_dem, canonical
from simulation_methods import build_mixed_stan_data, to_msprime_demography

POPS = ["a", "b", "c"]
SCENARIOS = ["hidden_loop", "no_loop"]
SOURCES = ["true", "hapibd"]
VARIANTS = {
    f"{likelihood}_{ne}": f"mixed_model_Nsmooth_{likelihood}"
    + ("_separate_ne" if ne == "separate" else "") + ".stan"
    for likelihood in ("poisson", "normal") for ne in ("shared", "separate")
}
UPPER = np.triu_indices(3)


def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def load_config(path=HERE / "config.json"):
    c = json.loads(Path(path).read_text())
    lo, hi, step = c["bin_edges_cm"]
    assert lo > 0 and hi > lo and step > 0
    assert np.isclose((hi-lo)/step, round((hi-lo)/step))
    assert c["blocks"] >= 2 and c["replicates"] >= 1 and c["pools"] >= 1
    assert c["block_cm"] > 0 and c["genome_cm"] > 0
    n = c["genome_cm"] / c["block_cm"]
    assert np.isclose(n, round(n)) and 2 <= round(n) <= c["blocks"]
    assert c["diploid_samples"] >= 2 and c["haploid_ne"] > 0
    assert c["recombination_rate"] > 0 and c["mutation_rate"] > 0
    assert c["ancestry_model"] in ("smc", "smc_prime", "hudson")
    assert 0 <= c["snp_min_maf"] < 0.5
    assert c["snp_se_block_cm"] > 0
    assert np.isclose(c["block_cm"] / c["snp_se_block_cm"],
                      round(c["block_cm"] / c["snp_se_block_cm"]))
    times = [c["times"][k] for k in ("c_split", "b_split", "left", "c_close", "right", "root")]
    assert all(b-a > 1 for a, b in zip([0] + times[:-1], times))
    assert c["T_max"] > times[-1]
    assert all(0 < f < 1 for f in c["fractions"].values())
    assert c["hapibd_min_output"] <= lo
    assert c["fit_seeds"] and c["pathfinder_paths"] > 0
    assert c["pathfinder_draws"] >= c["pathfinder_paths"]
    return c


def edges(c):
    lo, hi, step = c["bin_edges_cm"]
    return np.linspace(lo, hi, round((hi-lo)/step) + 1)


def pair_counts(c):
    n = 2 * c["diploid_samples"]
    return np.full((3, 3), n*n) - np.eye(3, dtype=int) * (n*n - n*(n-1)//2)


def simulation_demography(c, scenario):
    d = DemographicTopology(POPS)
    t = c["times"]
    if scenario == "hidden_loop":
        d.add_admixture_event("c", "c1", "c2")
    elif scenario != "no_loop":
        raise ValueError(scenario)
    d.add_admixture_event("b", "b1", "b2")
    d.add_merge_event("a", "b1", "left")
    if scenario == "hidden_loop":
        d.add_merge_event("c1", "c2", "anc_c")
    d.add_merge_event("b2", "anc_c" if scenario == "hidden_loop" else "c", "right")
    d.add_merge_event("left", "right", "root")
    for name in d.nodes:
        d.set_node_ne(name, c["haploid_ne"] / 2)
    if scenario == "hidden_loop":
        d.set_admixture_parameters("c", t["c_split"], c["fractions"]["c"], "c1")
        d.set_merge_time("anc_c", t["c_close"])
    d.set_admixture_parameters("b", t["b_split"], c["fractions"]["b"], "b1")
    for name in ("left", "right", "root"):
        d.set_merge_time(name, t[name])
    d.finalize_root()
    return d


def replay(events):
    d = DemographicTopology(POPS)
    for ev in events:
        if ev["type"] == "MERGE":
            d.add_merge_event(*ev["children"], ev["parent"])
        else:
            d.add_admixture_event(ev["child"], *ev["parents"])
    return d


def candidates():
    """21 shapes, with ALL valid chronological orders, including pre-admix merges."""
    result = []
    for graph in enumerate_all(POPS):
        base = quiet(make_dem, POPS, graph)
        valid = []
        for order in itertools.permutations(base.ordered_events):
            alive = set(POPS)
            for ev in order:
                children = set(ev["children"] if ev["type"] == "MERGE" else [ev["child"]])
                if not children <= alive:
                    break
                alive -= children
                alive.update([ev["parent"]] if ev["type"] == "MERGE" else ev["parents"])
            else:
                if len(alive) == 1:
                    valid.append(order)
        for j, order in enumerate(valid):
            d = quiet(replay, order)
            result.append({"id": f"g{graph['index']:02d}_o{j+1}",
                           "graph": graph["index"], "newick": graph["newick"],
                           "admixed": graph["admixed"], "n_orders": len(valid),
                           "events": d.ordered_events, "nodes": list(d.nodes),
                           "correct": graph["admixed"] == "b" and graph["newick"] ==
                           canonical((("a", "b.1"), ("b.2", "c")))})
    return result


def smooth_data(d):
    idx = {name: i for i, name in enumerate(d.nodes)}
    parent, admix, start = ([0]*len(idx) for _ in range(3))
    a = 0
    for k, ev in enumerate(d.ordered_events, 1):
        if ev["type"] == "MERGE":
            start[idx[ev["parent"]]] = k
            for child in ev["children"]:
                parent[idx[child]] = idx[ev["parent"]] + 1
        else:
            a += 1
            admix[idx[ev["child"]]] = a
            for name in ev["parents"]:
                start[idx[name]] = k
    return dict(ne_parent=parent, ne_admix_idx=admix, ne_start_event=start)


def selections(c, pool):
    # Same block IDs across sources/models/scenarios; scenarios have different ancestry seeds.
    rng = np.random.default_rng(np.random.SeedSequence([c["seed"], pool, 987]))
    n = round(c["genome_cm"] / c["block_cm"])
    return [sorted(rng.choice(c["blocks"], n, replace=False).tolist())
            for _ in range(c["replicates"])]


def sample_pairs(ts):
    pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
    lookup = {int(u): i for i, name in enumerate(POPS) for u in ts.samples(population=pop_id[name])}
    uv = np.asarray(list(itertools.combinations(ts.samples(), 2)), dtype=np.int32)
    cells = np.asarray([min(lookup[int(u)], lookup[int(v)])*3 + max(lookup[int(u)], lookup[int(v)])
                        for u, v in uv], dtype=np.int32)
    return uv, cells, lookup


def mrca_library():
    folder = HERE / ".build"
    folder.mkdir(exist_ok=True)
    src = HERE / "mrca_scan.cpp"
    lib = folder / ("mrca_" + hashlib.sha256(src.read_bytes()).hexdigest()[:12] + ".so")
    if not lib.exists():
        compiler = shutil.which(os.environ.get("CXX", "c++"))
        if not compiler:
            raise RuntimeError("C++ compiler required for fast strict-MRCA extraction (CXX or c++).")
        tmp = lib.with_suffix(".tmp.so")
        subprocess.run([compiler, "-O3", "-std=c++11", "-shared", "-fPIC", str(src), "-o", str(tmp)], check=True)
        tmp.replace(lib)
    fn = ctypes.CDLL(str(lib)).scan
    ip = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
    dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    lp = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
    fn.argtypes = [ctypes.c_int, ip, ip, ip, ip, dp, ctypes.c_double,
                   ctypes.c_double, ctypes.c_int, dp, ip, dp, lp, dp, ctypes.c_int]
    fn.restype = None
    return fn


def true_ibd(ts, c):
    uv, cells, _ = sample_pairs(ts)
    u, v = (np.ascontiguousarray(uv[:, k]) for k in (0, 1))
    previous = np.full(len(uv), -2, np.int32)
    start = np.zeros(len(uv))
    e = edges(c)
    count = np.zeros((len(e)-1, 3, 3), np.int64)
    length = np.zeros_like(count, dtype=float)
    fn = mrca_library()
    times = np.ascontiguousarray(ts.nodes_time)
    for tree in ts.trees():
        fn(len(u), u, v, cells, tree.parent_array, times, tree.interval.left,
           100*c["recombination_rate"], len(e)-1, e, previous, start, count, length, 0)
    fn(len(u), u, v, cells, tree.parent_array, times, ts.sequence_length,
       100*c["recombination_rate"], len(e)-1, e, previous, start, count, length, 1)
    for i, j in zip(*UPPER):
        count[:, j, i], length[:, j, i] = count[:, i, j], length[:, i, j]
    return count, length


def hapibd(ts, c, folder, command):
    """command is an argv prefix, e.g. ['java','-jar','/path/hap-ibd.jar']."""
    import gzip
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    vcf, gmap, out = folder / "input.vcf", folder / "input.map", folder / "hap"
    individuals = [ind for ind in ts.individuals() if len(ind.nodes) == 2 and
                   all(ts.node(int(u)).is_sample() for u in ind.nodes)]
    names = [f"sample_{ind.id}" for ind in individuals]
    mapping = {(name, k+1): int(u) for name, ind in zip(names, individuals)
               for k, u in enumerate(ind.nodes)}
    if len(mapping) != ts.num_samples:
        raise ValueError("Every sample must belong to a diploid individual.")
    with vcf.open("w") as handle:
        ts.write_vcf(handle, contig_id="1", individuals=[x.id for x in individuals],
                     individual_names=names, position_transform=lambda x: np.asarray(x)+1)
    scale = 100*c["recombination_rate"]
    # +1 VCF positions, exact same slope as the simulation (cM/bp).
    gmap.write_text(f"1 . 0 1\n1 . {ts.sequence_length*scale:.12g} {int(ts.sequence_length)+1}\n")
    args = list(command) + [f"gt={vcf}", f"map={gmap}", f"out={out}",
                            f"min-seed={c['hapibd_min_seed']}",
                            f"min-output={c['hapibd_min_output']}",
                            f"nthreads={c['hapibd_threads']}"]
    with (folder / "command.log").open("w") as handle:
        handle.write(json.dumps(args) + "\n")
        handle.flush()
        subprocess.run(args, stdout=handle, stderr=subprocess.STDOUT, check=True)
    _, _, lookup = sample_pairs(ts)
    count = np.zeros((len(edges(c))-1, 3, 3), np.int64)
    length = np.zeros_like(count, dtype=float)
    # .hbd.gz contains the homologous pair within each diploid; the Stan
    # n*(n-1)/2 exposure includes these pairs, so BOTH files are required.
    for suffix in (".ibd.gz", ".hbd.gz"):
        with gzip.open(str(out)+suffix, "rt") as handle:
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
                b = np.searchsorted(edges(c), size, side="left") - 1
                if 0 <= b < len(count):
                    count[b, i, j] += 1
                    length[b, i, j] += size
    for i, j in zip(*UPPER):
        count[:, j, i], length[:, j, i] = count[:, i, j], length[:, i, j]
    return count, length


def snp_summaries(ts, c):
    """TreeMix ratio-of-sums, same ancestral ascertainment as old simulations.

    Store physical sub-block sufficient statistics rather than every SNP.
    Unlike the old SNP-count windows, SE blocks cannot bridge chromosomes.
    """
    pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
    samples = list(ts.samples())
    column = {int(u): j for j, u in enumerate(samples)}
    indices = [[column[int(u)] for u in ts.samples(population=pop_id[p])] for p in POPS]
    n = round(c["block_cm"] / c["snp_se_block_cm"])
    numer = np.zeros((n, 3, 3)); denom = np.zeros(n); sites = np.zeros(n, np.int64)
    for var in ts.variants():
        # Explicitly discard multiallelic/missing sites: genotype codes are not allele counts.
        if len(var.alleles) != 2 or np.any(var.genotypes < 0):
            continue
        if c["snp_ancestral_only"] and (len(var.site.mutations) != 1 or
                ts.node(var.site.mutations[0].node).time < c["times"]["root"]):
            continue
        freq = np.asarray([np.mean(var.genotypes[ix]) for ix in indices])
        mean = freq.mean()
        if not c["snp_min_maf"] < mean < 1-c["snp_min_maf"]:
            continue
        k = min(n-1, int(var.site.position * 100*c["recombination_rate"] / c["snp_se_block_cm"]))
        dev = freq - mean
        numer[k] += np.outer(dev, dev)
        denom[k] += mean*(1-mean)
        sites[k] += 1
    return dict(snp_numer=numer, snp_denom=denom, snp_sites=sites)


def snp_covariance(numer, denom, c):
    numer = np.asarray(numer).reshape(-1, 3, 3)
    denom = np.asarray(denom).reshape(-1)
    if len(denom) < 2 or np.sum(denom) <= 0 or np.any(denom.sum()-denom <= 0):
        raise ValueError("Not enough SNP-bearing blocks to estimate uncertainty")
    # Centered finite-sample correction retained from inference_methods.py.
    center = np.eye(3) - np.ones((3, 3))/3
    correction = center @ np.diag(np.full(3, 1/(2*c["diploid_samples"]))) @ center
    w = numer.sum(axis=0) / denom.sum() - correction
    leave = (numer.sum(axis=0)-numer) / (denom.sum()-denom)[:, None, None] - correction
    se = np.sqrt((len(denom)-1)/len(denom)*np.sum((leave-leave.mean(axis=0))**2, axis=0))
    return w, np.maximum(se, 1e-8)


def aggregate(blocks, chosen, source, c):
    count = np.sum([blocks[i][source+"_count"] for i in chosen], axis=0)
    lengths = np.asarray([blocks[i][source+"_length"] for i in chosen])
    cm = len(chosen)*c["block_cm"]
    fraction = lengths.sum(axis=0) / (cm * pair_counts(c))
    per_block = lengths / (c["block_cm"] * pair_counts(c))
    se = np.maximum(per_block.std(axis=0, ddof=1)/np.sqrt(len(chosen)), 1e-8)
    w, ws = snp_covariance([blocks[i]["snp_numer"] for i in chosen],
                           [blocks[i]["snp_denom"] for i in chosen], c)
    return dict(ibd_count=count, ibd_hat=fraction, ibd_se=se, w_hat=w, w_se=ws, cm=cm)


def stan_data(candidate, obs, c):
    d = quiet(replay, candidate["events"])
    e = edges(c); bins = list(zip(e[:-1], e[1:]))
    data = quiet(build_mixed_stan_data, d,
                 dict(enumerate(obs["ibd_hat"])), dict(enumerate(obs["ibd_se"]**2)), bins,
                 obs["w_hat"], obs["w_se"], n_samples_per_pop=2*c["diploid_samples"],
                 T_max=c["T_max"], cm=obs["cm"], ibd_count=dict(enumerate(obs["ibd_count"])))
    # Empty arrays must retain their intended 0 x 4 / 0 x n_nodes x n_nodes shape.
    data["admixture_map"] = np.asarray(data["admixture_map"], dtype=int).reshape(-1, 4)
    data.update(smooth_data(d))
    return data
