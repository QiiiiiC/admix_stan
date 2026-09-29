"""b-loop then admixture, with exponential growth on a and b and 3:1 branch-size
asymmetries in the loop and in the admixture.  Data, generating history and
candidate graphs.

Every fitted candidate has a CONSTANT Ne per branch (Nsmooth), so under the
`growth` scenario every candidate -- including the explicit-loop graph 22 --
is misspecified in Ne.  The `constant` control has the identical topology and
times with all branches at haploid Ne 15,000, so it isolates what the Ne
misspecification alone does to topology selection and to shared-vs-separate Ne.

Generic machinery (simulation, IBD/SNP extraction, block aggregation, Stan data,
models, fitting) lives in the repository's `methods` package; this file holds
only what is specific to this study.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".build" / "matplotlib"))

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from methods import blocks, ibd, snp
from methods import stan_data as stan
from methods.demography import DemographicTopology, replay as _replay, clades as _clades
from methods.diagnostics import residual_tilt as _residual_tilt, trajectory_reference
from methods.models import MODELS
from methods.simulate import add_branch, size_at as _size_at
from methods.stan_data import smooth_data  # noqa: F401  (re-exported for tests)
from methods.topologies import canonical, ordered_candidates
from methods.utils import quiet, save_json, digest  # noqa: F401  (re-exported for the runner)

POPS = ["a", "b", "c"]
SCENARIOS = ["growth", "constant"]
SOURCES = ["true", "hapibd"]
VARIANTS = {"poisson_shared": MODELS["mixed_Nsmooth_poisson"],
            "poisson_separate": MODELS["mixed_Nsmooth_poisson_separate"]}
UPPER = np.triu_indices(3)

# Generating branches as msprime populations (no dots: msprime names must be identifiers).
GENERATING = ["a", "b", "c", "loop1", "loop2", "anc_b", "b1", "b2", "left", "right", "root"]
# Branches a loop-omitting candidate can name.  Its "b" spans the leaf, the loop
# and anc_b, so its true size is a COMPOSITE -- see effective_ne.
BACKBONE = ["a", "b", "c", "b.1", "b.2", "left", "right", "root"]


def load_config(path=HERE / "config.json"):
    import json
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
    return ibd.bin_edges(*c["bin_edges_cm"])


def n_haploid(c):
    return [2 * c["diploid_samples"]] * len(POPS)


def pair_counts(c):
    return ibd.pair_counts(n_haploid(c))


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
    return _size_at(*branch_sizes(c, scenario)[name], *branch_interval(c, name), t)


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
    if name == "b" and collapsed:
        t0, t1 = 0, c["times"]["b_split"]
    else:
        t0, t1 = branch_interval(c, {"b.1": "b1", "b.2": "b2"}.get(name, name))
    ref = trajectory_reference(lambda t: effective_ne(c, scenario, name, t, collapsed), t0, t1)
    return dict(branch=name, **ref)


def msprime_demography(c, scenario):
    """The generating history in msprime, with per-branch growth (sizes haploid,
    halved inside `add_branch`).  Verified against DemographyDebugger in tests."""
    import msprime
    d = msprime.Demography()
    sizes = branch_sizes(c, scenario)
    for name in GENERATING:
        add_branch(d, name, *sizes[name], *branch_interval(c, name))
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
    return _replay(POPS, events)


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
    backbone = canonical((("a", "b.1"), ("b.2", "c")))
    result = [dict(x, explicit_loop=False, correct=x["admixed"] == "b" and x["newick"] == backbone)
              for x in ordered_candidates(POPS)]
    for j, merge_order in enumerate((("left", "right"), ("right", "left")), 1):
        d = _explicit_loop(merge_order)
        result.append(dict(id=f"g22_o{j}", graph=22, newick="b-loop + ((a,b.1),(b.2,c))",
                           admixed="b", n_orders=2, events=d.ordered_events, nodes=list(d.nodes),
                           correct=True, explicit_loop=True))
    return result


def clades(candidate):
    return _clades(POPS, candidate["events"])


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


def selections(c, pool):
    # Same block IDs across sources/models/scenarios; scenarios have different ancestry seeds.
    return blocks.selections(c["seed"], pool, c["blocks"], round(c["genome_cm"] / c["block_cm"]),
                             c["replicates"])


# ---------------------------------------------------------------------------
# Observations, with this study's configuration bound in
# ---------------------------------------------------------------------------
def true_ibd(ts, c):
    return ibd.true_ibd(ts, POPS, edges(c), c["recombination_rate"])


def hapibd(ts, c, folder, command):
    """command is an argv prefix, e.g. ['java','-jar','/path/hap-ibd.jar']."""
    return ibd.hapibd(ts, POPS, edges(c), c["recombination_rate"], folder, command,
                      c["hapibd_min_seed"], c["hapibd_min_output"], c["hapibd_threads"])


def snp_summaries(ts, c):
    """TreeMix ratio-of-sums with ancestral-SNP ascertainment; per-sub-block sums."""
    return snp.snp_summaries(ts, POPS, c["block_cm"], c["snp_se_block_cm"], c["recombination_rate"],
                             c["snp_min_maf"], c["times"]["root"] if c["snp_ancestral_only"] else None)


def aggregate(block_list, chosen, source, c):
    return blocks.aggregate(block_list, chosen, source, c["block_cm"], n_haploid(c))


def stan_data(candidate, obs, c, variant="poisson_shared"):
    return stan.stan_data(replay(candidate["events"]), obs, VARIANTS[variant],
                          edges=edges(c), n_haploid=n_haploid(c), T_max=c["T_max"])


def residual_tilt(pearson, c):
    return _residual_tilt(pearson, edges(c))
