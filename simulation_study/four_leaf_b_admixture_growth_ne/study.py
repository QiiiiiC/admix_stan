"""Four leaves, one admixture, growth-like branch sizes: topology, recovery and ELBO vs log Z.

Generating history (backwards in time):

  20   b splits into b1 (0.7) and b2 (0.3)         admixture
  60   a and b1 merge (-> ab)
  100  b2 and c merge (-> cb)
  200  cb and d merge (-> cbd)
  400  ab and cbd merge (-> root)

Branch sizes come from the config, either as a rule or as a table (`branch_ne`):
  "growth"     mimic exponential growth while staying constant along each branch:
               every lineage grows log-linearly from the ancestral size at the root to
               its own present-day size; a branch shared by several leaves follows the
               geometric mean of their present sizes; each branch takes that curve's
               value at its time midpoint.  Recent branches are large and different.
  "haploid_ne" an explicit haploid Ne per branch (the four_leaf_b_admixture_varying_ne
               design, rerun through this pipeline from its saved tree sequences).

Candidates (a graph is an event order): T_true, and T_null = ((a,b),(c,d)) with b
wholly on its major source side.  Models, both Nsmooth with the Poisson IBD likelihood:
`mixed` (IBD + SNP) and `ibd` (IBD only).  IBD sources (config `sources`): `true`
(segments read off the simulated genealogies) and `hapibd` (hap-IBD called on the
simulated genotypes).  NUTS runs on one source (config `nuts_source`).

Generic machinery lives in the repository's `methods` package.
"""
from __future__ import annotations

import json
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
from methods.demography import DemographicTopology, replay as _replay
from methods.models import MODELS
from methods.simulate import to_msprime
from methods.utils import quiet, save_json, digest  # noqa: F401  (re-exported for the runner)

POPS = ["a", "b", "c", "d"]
SOURCES = ["true", "hapibd"]
MODEL_NAMES = {"mixed": MODELS["mixed_Nsmooth_poisson"], "ibd": MODELS["ibd_Nsmooth_poisson"]}
NODES = ["a", "b", "c", "d", "b1", "b2", "ab", "cb", "cbd", "root"]
EVENTS = ["b_admix", "ab", "cb", "cbd", "root"]
# Leaves below each branch (admixture sources belong to b's lineage)
CLADE = {"a": "a", "b": "b", "c": "c", "d": "d", "b1": "b", "b2": "b", "ab": "ab", "cb": "bc", "cbd": "bcd"}


def load_config(path=HERE / "config.json"):
    c = json.loads(Path(path).read_text())
    e = np.asarray(c["bin_edges_cm"], float)
    assert len(e) >= 2 and np.all(np.diff(e) > 0) and e[0] > 0 and e[-1] <= c["block_cm"]
    for cm in c["genome_cm"]:
        n = cm / c["block_cm"]
        assert np.isclose(n, round(n)) and 2 <= round(n) <= c["blocks"], f"{cm} cM: whole blocks, 2..pool"
    assert c["genome_cm"] == sorted(set(c["genome_cm"])) and c["replicates"] >= c["nuts_replicates"] >= 1
    t = c["times"]
    assert list(t) == EVENTS and all(t[x] < t[y] for x, y in zip(EVENTS, EVENTS[1:])) and t["b_admix"] > 1
    assert c["T_max"] > t["root"] and 0 < c["fraction_b1"] < 1
    assert ("growth" in c) != ("haploid_ne" in c), "exactly one of growth / haploid_ne"
    if "growth" in c:
        g = c["growth"]
        assert sorted(g["present_ne"]) == POPS and all(v > g["ancestral_ne"] > 0 for v in g["present_ne"].values())
    else:
        assert sorted(c["haploid_ne"]) == sorted(NODES) and all(v > 0 for v in c["haploid_ne"].values())
    assert c["sources"] and set(c["sources"]) <= set(SOURCES) and c["nuts_source"] in c["sources"]
    if "hapibd" in c["sources"]:
        h = c["hapibd"]
        assert 0 < h["min_seed"] and 0 < h["min_output"] <= e[0], "hap-IBD min-output must not truncate the lowest bin"
    if "trees_from" in c:
        assert (ROOT / c["trees_from"]).is_dir(), "trees_from must name an existing run pool (relative to the repository root)"
    assert np.isclose(c["block_cm"] / c["snp_se_block_cm"], round(c["block_cm"] / c["snp_se_block_cm"]))
    return c


def edges(c):
    return np.asarray(c["bin_edges_cm"], float)


def n_haploid(c):
    return [2 * c["diploid_samples"]] * len(POPS)


# ---------------------------------------------------------------------------
# Generating history
# ---------------------------------------------------------------------------
def interval(c, node):
    t = c["times"]
    return {"a": (0, t["ab"]), "b": (0, t["b_admix"]), "c": (0, t["cb"]), "d": (0, t["cbd"]),
            "b1": (t["b_admix"], t["ab"]), "b2": (t["b_admix"], t["cb"]), "ab": (t["ab"], t["root"]),
            "cb": (t["cb"], t["cbd"]), "cbd": (t["cbd"], t["root"]), "root": (t["root"], np.inf)}[node]


def branch_ne(c):
    """Haploid Ne per branch.  Table designs return the table.  Growth designs: the
    log-linear growth curve of the branch's clade, from the ancestral size at the root
    to the clade's present size (geometric mean of its leaves), at the branch's time
    midpoint."""
    if "haploid_ne" in c:
        return {node: float(c["haploid_ne"][node]) for node in NODES}
    g = c["growth"]; t_root = c["times"]["root"]; n0 = g["ancestral_ne"]
    out = {}
    for node in NODES:
        if node == "root":
            out[node] = float(n0); continue
        t0, t1 = interval(c, node)
        present = np.exp(np.mean([np.log(g["present_ne"][leaf]) for leaf in CLADE[node]]))
        s = (t0 + t1) / 2 / t_root
        out[node] = float(np.exp(s * np.log(n0) + (1 - s) * np.log(present)))
    return out


def build():
    d = DemographicTopology(POPS)
    d.add_admixture_event("b", "b1", "b2")
    d.add_merge_event("a", "b1", "ab")
    d.add_merge_event("b2", "c", "cb")
    d.add_merge_event("cb", "d", "cbd")
    d.add_merge_event("ab", "cbd", "root")
    return d


def build_null():
    """No admixture: b wholly on its major (70%) source side, ((a,b),(c,d)), events in
    the true time order (ab 60, cd 200, root 400)."""
    d = DemographicTopology(POPS)
    d.add_merge_event("a", "b", "ab")
    d.add_merge_event("c", "d", "cd")
    d.add_merge_event("ab", "cd", "root")
    return d


GRAPHS = {"T_true": build, "T_null": build_null}


def true_topology(c):
    d = build()
    t = c["times"]
    for node, ne in branch_ne(c).items():
        d.set_node_ne(node, round(ne) // 2)          # the DSL stores diploid sizes
    d.set_admixture_parameters("b", t["b_admix"], c["fraction_b1"], "b1")
    for node in ("ab", "cb", "cbd", "root"):
        d.set_merge_time(node, t[node])
    d.finalize_root()
    d.is_valid()
    return d


def msprime_demography(c):
    return to_msprime(true_topology(c))


def parameters(c):
    """All 16 true parameters in Stan's indexing (T_true)."""
    truth = true_topology(c)
    out = [dict(name=f"t_{e}", variable="cumulative_times", index=k, truth=c["times"][e]) for k, e in enumerate(EVENTS)]
    out.append(dict(name="f_b1", variable="admixture_fractions", index=0, truth=c["fraction_b1"]))
    out += [dict(name=f"Ne_{n}", variable="Ne", index=k, truth=2 * truth.nodes[n].ne) for k, n in enumerate(NODES)]
    return out


# ---------------------------------------------------------------------------
# Subsamples, tasks, observations
# ---------------------------------------------------------------------------
def selections(c, pool, cm):
    return blocks.selections(c["seed"], pool, c["blocks"], round(cm / c["block_cm"]), c["replicates"], int(cm))


def fit_tasks(c, sources, models, graphs, replicate_limit=None):
    reps = range(min(c["replicates"], replicate_limit or c["replicates"]))
    return [(cm, rep, src, m, g) for cm in c["genome_cm"] for rep in reps for src in sources for m in models
            for g in graphs]


def nuts_tasks(c, graphs):
    return [(cm, rep, c["nuts_source"], "mixed", g) for cm in c["genome_cm"] for rep in range(c["nuts_replicates"])
            for g in graphs]


def true_ibd(ts, c):
    return ibd.true_ibd(ts, POPS, edges(c), c["recombination_rate"])


def hapibd(ts, c, folder, command):
    h = c["hapibd"]
    return ibd.hapibd(ts, POPS, edges(c), c["recombination_rate"], folder, command, h["min_seed"], h["min_output"],
                      h["threads"])


def snp_summaries(ts, c):
    return snp.snp_summaries(ts, POPS, c["block_cm"], c["snp_se_block_cm"], c["recombination_rate"],
                             c["snp_min_maf"], c["times"]["root"] if c["snp_ancestral_only"] else None)


def aggregate(block_list, chosen, c, source):
    return blocks.aggregate(block_list, chosen, source, c["block_cm"], n_haploid(c))


def stan_data(graph, obs, c, model):
    return stan.stan_data(GRAPHS[graph](), obs, MODEL_NAMES[model], edges=edges(c), n_haploid=n_haploid(c),
                          T_max=c["T_max"])


def replay(events):
    return _replay(POPS, events)
