"""b-loop then admixture, with exponential growth on a and b and 3:1 branch-size
asymmetries in the loop and in the admixture.  Data, generating history and
candidate graphs.

Every fitted candidate has a CONSTANT Ne per branch (Nsmooth), so under the
`growth` scenario every candidate -- including the explicit-loop graph 22 --
is misspecified in Ne.  The `constant` control has the identical topology and
times with all branches at haploid Ne 15,000, so it isolates what the Ne
misspecification alone does to topology selection and to shared-vs-separate Ne.
"""
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
from simulation_methods import build_mixed_stan_data

POPS = ["a", "b", "c"]
SCENARIOS = ["growth", "constant"]
SOURCES = ["true", "hapibd"]
VARIANTS = {
    f"{likelihood}_{ne}": f"mixed_model_Nsmooth_{likelihood}"
    + ("_separate_ne" if ne == "separate" else "") + ".stan"
    for likelihood in ("poisson",) for ne in ("shared", "separate")
}
UPPER = np.triu_indices(3)

# Generating branches as msprime populations (no dots: msprime names must be identifiers).
GENERATING = ["a", "b", "c", "loop1", "loop2", "anc_b", "b1", "b2", "left", "right", "root"]
# Branches a loop-omitting candidate can name.  Its "b" spans the leaf, the loop
# and anc_b, so its true size is a COMPOSITE -- see effective_ne.
BACKBONE = ["a", "b", "c", "b.1", "b.2", "left", "right", "root"]


def quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def save_json(path, value):
    """Atomic write.  The temp name carries the PID: parallel workers write
    identical provenance/manifest files, and a shared temp path let one worker's
    rename delete another's file mid-flight (FileNotFoundError on replace)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f"{path.suffix}.{os.getpid()}.tmp")
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
    times = [c["times"][k] for k in ("loop_open", "loop_close", "b_split", "left", "right", "root")]
    assert all(b-a > 1 for a, b in zip([0] + times[:-1], times))
    assert c["T_max"] > times[-1]
    assert all(0 < f < 1 for f in c["fractions"].values())
    assert all(c["ne"][k] >= 1 for k in ("growth_fold_a", "growth_fold_b", "ratio_loop", "ratio_split"))
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


# ---------------------------------------------------------------------------
# Generating history: intervals, sizes, growth
# ---------------------------------------------------------------------------
def branch_interval(c, name):
    """(young end, old end) in generations for a generating branch."""
    t = c["times"]
    return {"a": (0, t["left"]), "b": (0, t["loop_open"]), "c": (0, t["right"]),
            "loop1": (t["loop_open"], t["loop_close"]), "loop2": (t["loop_open"], t["loop_close"]),
            "anc_b": (t["loop_close"], t["b_split"]),
            "b1": (t["b_split"], t["left"]), "b2": (t["b_split"], t["right"]),
            "left": (t["left"], t["root"]), "right": (t["right"], t["root"]),
            "root": (t["root"], np.inf)}[name]


def branch_sizes(c, scenario):
    """Haploid Ne at the (young, old) ends of every generating branch.

    growth  : a grows growth_fold_a-fold toward the present across its whole leaf
              branch; b grows growth_fold_b-fold across its leaf branch (present
              back to loop opening); loop1 = ratio_loop x loop2; b1 = ratio_split
              x b2.  The larger side of each ratio sits at haploid_ne, so the
              minor side is the small one.  All other branches are haploid_ne.
    constant: every branch haploid_ne.  Same topology, same times.
    """
    n0 = c["haploid_ne"]; g = c["ne"]
    sizes = {name: (n0, n0) for name in GENERATING}
    if scenario == "growth":
        sizes["a"] = (n0 * g["growth_fold_a"], n0)
        sizes["b"] = (n0 * g["growth_fold_b"], n0)
        sizes["loop2"] = (n0 / g["ratio_loop"], n0 / g["ratio_loop"])
        sizes["b2"] = (n0 / g["ratio_split"], n0 / g["ratio_split"])
    elif scenario != "constant":
        raise ValueError(scenario)
    return sizes


def size_at(c, scenario, name, t):
    """Haploid size of one generating branch at time t (exponential between ends)."""
    young, old = branch_sizes(c, scenario)[name]
    t0, t1 = branch_interval(c, name)
    t = np.asarray(t, float)
    if young == old or not np.isfinite(t1):
        return np.full(t.shape, float(young))
    return young * np.exp(np.log(old/young) * np.clip((t-t0)/(t1-t0), 0, 1))


def effective_ne(c, scenario, name, t, collapsed=True):
    """Pairwise-coalescence effective size of a branch a candidate can name.

    For the collapsed backbone's "b" (present back to the external admixture)
    the loop interval holds two branches; two b lineages fall into branch k
    independently with probability p_k, so the coalescence rate is
    sum_k p_k^2 / N_k and the effective size is its reciprocal.  The same
    quantity governs allele-frequency drift of the admixed population
    (variance = sum_k p_k^2 var_k), so one reference serves both data types.
    """
    t = np.asarray(t, float)
    if name == "b" and collapsed:
        f = c["fractions"]["loop"]; tt = c["times"]
        leaf = size_at(c, scenario, "b", t)
        loop = 1.0 / (f**2 / size_at(c, scenario, "loop1", t) + (1-f)**2 / size_at(c, scenario, "loop2", t))
        anc = size_at(c, scenario, "anc_b", t)
        return np.where(t < tt["loop_open"], leaf, np.where(t < tt["loop_close"], loop, anc))
    return size_at(c, scenario, {"b.1": "b1", "b.2": "b2"}.get(name, name), t)


def reference(c, scenario, name, collapsed=True):
    """What a constant-Ne branch could hope to recover: the ends, the harmonic
    mean (total drift over the branch is duration / harmonic mean, so this is
    the SNP-side target) and the arithmetic mean, over the true interval."""
    t = c["times"]
    if name == "b" and collapsed:
        t0, t1 = 0, t["b_split"]
    else:
        t0, t1 = branch_interval(c, {"b.1": "b1", "b.2": "b2"}.get(name, name))
    if not np.isfinite(t1):
        n = float(effective_ne(c, scenario, name, t0, collapsed))
        return dict(branch=name, t0=t0, t1=None, young=n, old=n, harmonic=n, arithmetic=n)
    grid = np.linspace(t0, t1, 2001)
    n = effective_ne(c, scenario, name, grid, collapsed)
    return dict(branch=name, t0=t0, t1=t1, young=float(n[0]), old=float(n[-1]),
                harmonic=float((t1-t0) / np.trapz(1/n, grid)),
                arithmetic=float(np.trapz(n, grid) / (t1-t0)))


def msprime_demography(c, scenario):
    """Build the generating history directly in msprime, with per-branch growth.

    msprime measures growth from time 0, so an ancestral branch that only exists
    on [t0, t1] gets a population_parameters_change AT t0 setting its young-end
    size and rate; that anchors size(t0) exactly.  Growth is switched off at t1
    so no branch keeps shrinking after its lineages have moved on.  Sizes are
    passed as diploid (haploid / 2).  Verified against DemographyDebugger.
    """
    import msprime
    d = msprime.Demography()
    sizes = branch_sizes(c, scenario)
    for name in GENERATING:
        young, old = sizes[name]
        t0, t1 = branch_interval(c, name)
        rate = 0.0 if (young == old or not np.isfinite(t1)) else float(np.log(young/old) / (t1-t0))
        d.add_population(name=name, initial_size=young/2, growth_rate=rate if t0 == 0 else 0.0)
        if rate and t0 > 0:
            d.add_population_parameters_change(time=t0, population=name, initial_size=young/2, growth_rate=rate)
        if rate:
            d.add_population_parameters_change(time=t1, population=name, growth_rate=0.0)
    t = c["times"]; f = c["fractions"]
    d.add_admixture(time=t["loop_open"], derived="b", ancestral=["loop1", "loop2"],
                    proportions=[f["loop"], 1-f["loop"]])
    d.add_population_split(time=t["loop_close"], derived=["loop1", "loop2"], ancestral="anc_b")
    d.add_admixture(time=t["b_split"], derived="anc_b", ancestral=["b1", "b2"],
                    proportions=[f["b"], 1-f["b"]])
    d.add_population_split(time=t["left"], derived=["a", "b1"], ancestral="left")
    d.add_population_split(time=t["right"], derived=["b2", "c"], ancestral="right")
    d.add_population_split(time=t["root"], derived=["left", "right"], ancestral="root")
    d.sort_events()
    d.validate()
    return d


# ---------------------------------------------------------------------------
# Candidate graphs (identical search to b_loop_then_admixture)
# ---------------------------------------------------------------------------
def replay(events):
    d = DemographicTopology(POPS)
    for ev in events:
        if ev["type"] == "MERGE":
            d.add_merge_event(*ev["children"], ev["parent"])
        else:
            d.add_admixture_event(ev["child"], *ev["parents"])
    return d


def _explicit_loop(merge_order):
    d = DemographicTopology(POPS)
    d.add_admixture_event("b", "loop1", "loop2")
    d.add_merge_event("loop1", "loop2", "anc_b")
    d.add_admixture_event("anc_b", "b.1", "b.2")
    for side in merge_order:
        if side == "left":
            d.add_merge_event("a", "b.1", "left")
        else:
            d.add_merge_event("b.2", "c", "right")
    d.add_merge_event("left", "right", "root")
    return d


def candidates():
    """21 omitted-loop shapes + explicit generating loop graph; 29 total orders."""
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
                           "explicit_loop": False,
                           "correct": graph["admixed"] == "b" and graph["newick"] ==
                           canonical((("a", "b.1"), ("b.2", "c")))})
    for j, merge_order in enumerate((("left", "right"), ("right", "left")), 1):
        d = quiet(_explicit_loop, merge_order)
        result.append(dict(id=f"g22_o{j}", graph=22, newick="b-loop + ((a,b.1),(b.2,c))",
                           admixed="b", n_orders=2, events=d.ordered_events, nodes=list(d.nodes),
                           correct=True, explicit_loop=True))
    return result


def clades(candidate):
    """Leaf-set under every node, walking the candidate's events."""
    out = {name: {name} for name in POPS}
    for ev in candidate["events"]:
        if ev["type"] == "ADMIXTURE":
            for p in ev["parents"]:
                out[p] = {p}
        else:
            out[ev["parent"]] = set.union(*(out[p] for p in ev["children"]))
    return out


def branch_map(candidate):
    """Candidate node -> generating/backbone branch it stands for, for correct
    candidates only.  Internal nodes of the omitted-loop backbone are named by
    what they contain.  Returns {} for wrong graphs (no comparable branch)."""
    if not candidate["correct"]:
        return {}
    cl = clades(candidate)
    out = {}
    for node in candidate["nodes"]:
        members = cl[node]
        if node in ("a", "b", "c", "b.1", "b.2", "loop1", "loop2", "anc_b"):
            out[node] = node
        elif members == {"loop1", "loop2"}:
            out[node] = "anc_b"
        elif {"a", "c"} <= members:
            out[node] = "root"
        elif "a" in members:
            out[node] = "left"
        else:
            out[node] = "right"
    return out


def event_parameters(candidate, c):
    """Semantic event labels, including the distinction between the TWO admixtures.

    For the explicit loop the external b admixture is fraction index 1; the
    loop fraction is index 0 and is source-exchange symmetric. Report its
    larger source fraction, rather than treating source labels as identifiable.
    """
    if not candidate["correct"]:
        return []
    cl = {name: {name} for name in POPS}
    mapping = []
    admixture_index = 0
    for k, ev in enumerate(candidate["events"]):
        if ev["type"] == "ADMIXTURE":
            is_loop = ev["parents"] == ["loop1", "loop2"]
            name = "loop_open" if is_loop else "b_split"
            mapping.append(dict(name=name, variable="cumulative_times", index=k,
                                truth=c["times"][name], fold=False))
            mapping.append(dict(name="loop_major_fraction" if is_loop else "b_fraction",
                                variable="admixture_fractions", index=admixture_index,
                                truth=max(c["fractions"]["loop"], 1-c["fractions"]["loop"]) if is_loop else c["fractions"]["b"],
                                fold=is_loop))
            admixture_index += 1
            for p in ev["parents"]: cl[p] = {p}
        else:
            members = set.union(*(cl[p] for p in ev["children"]))
            cl[ev["parent"]] = members
            name = ("loop_close" if members == {"loop1", "loop2"} else
                    "root" if {"a", "c"} <= members else "left" if "a" in members else "right")
            mapping.append(dict(name=name, variable="cumulative_times", index=k,
                                truth=c["times"][name], fold=False))
    return mapping


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


# ---------------------------------------------------------------------------
# Observations (unchanged from b_loop_then_admixture)
# ---------------------------------------------------------------------------
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
    # .hbd.gz holds the homologous pair inside each diploid; the Stan diagonal
    # exposure n*(n-1)/2 includes those pairs, so BOTH files are required.
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
    """TreeMix ratio-of-sums with ancestral-SNP ascertainment; per-sub-block sums."""
    pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
    samples = list(ts.samples())
    column = {int(u): j for j, u in enumerate(samples)}
    indices = [[column[int(u)] for u in ts.samples(population=pop_id[p])] for p in POPS]
    n = round(c["block_cm"] / c["snp_se_block_cm"])
    numer = np.zeros((n, 3, 3)); denom = np.zeros(n); sites = np.zeros(n, np.int64)
    for var in ts.variants():
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
    data["admixture_map"] = np.asarray(data["admixture_map"], dtype=int).reshape(-1, 4)
    data.update(smooth_data(d))
    return data


def residual_tilt(pearson, c):
    """Least-squares slope of the Pearson residual against segment length, per
    pair, in residual units per cM.  A monotone tilt is the fingerprint of a
    branch whose Ne changed within the fitted interval: a constant Ne cannot
    make both the short (older) and long (recent) bins agree at once."""
    x = (edges(c)[:-1] + edges(c)[1:]) / 2
    x = x - x.mean()
    out = {}
    for i, j in zip(*UPPER):
        y = np.asarray(pearson)[:, i, j]
        out[f"{i},{j}"] = dict(slope=float(np.sum(x*(y-y.mean())) / np.sum(x*x)),
                               correlation=float(np.corrcoef(x, y)[0, 1]) if np.std(y) > 0 else 0.0)
    return out
