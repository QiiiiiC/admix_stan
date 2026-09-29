"""Recent admix-through-b loop + ancient admixture, uniform Ne: the Nfixed
experiment in which the mixed model beats both single-data models.

Generating history (backwards in time, haploid Ne 10,000 on every branch):

  20   a splits into aP1 (0.5) and aP2 (0.5)          recent admixture on a
  60   aP2 merges into b (-> bP)                       donor branch through b
  100  bP and aP1 merge (-> ab)                        loop closes
  500  c splits into cP1 (0.7) and cP2 (0.3)          ancient admixture on c
  700  ab and cP1 merge (-> left)
  900  d and cP2 merge (-> right)
  1000 left and right merge (-> root)

Candidates, fitted with IBD-only, SNP-only and mixed Nfixed models (Poisson IBD
likelihood on segment counts):

  T_true  the generating graph
  T_alt1  no recent loop (a and b merge directly); ancient admixture kept
  T_alt2  recent loop kept; no ancient admixture (c joins (ab, d) as a tree)

The prediction is a division of labour: IBD sees the recent loop, SNP the
ancient admixture, and only the mixed model prefers T_true on both contrasts.
This study reruns the original (new_pipeline/Nfixed/
topology_admix_thru_b_plus_ancient) on the shared `methods` pipeline across
genome lengths.

Generic machinery lives in the repository's `methods` package; this file holds
only what is specific to this study.
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
from methods.diagnostics import residual_tilt  # noqa: F401  (re-exported)
from methods.models import MODELS
from methods.simulate import to_msprime
from methods.utils import quiet, save_json, digest  # noqa: F401  (re-exported for the runner)

POPS = ["a", "b", "c", "d"]
SOURCES = ["true", "hapibd"]
# Nfixed with the Poisson count likelihood for IBD (as in the latest studies): no
# IBD standard error enters, so short genomes with few blocks are well posed.
MODEL_NAMES = {"ibd": MODELS["ibd_Nfixed_poisson"], "snp": MODELS["snp_Nfixed"],
               "mixed": MODELS["mixed_Nfixed_poisson"]}
CANDIDATES = ["T_true", "T_alt1", "T_alt2"]
CONTRASTS = {"recent": ("T_true", "T_alt1"), "ancient": ("T_true", "T_alt2")}


def load_config(path=HERE / "config.json"):
    c = json.loads(Path(path).read_text())
    e = np.asarray(c["bin_edges_cm"], float)
    assert len(e) >= 2 and np.all(np.diff(e) > 0) and e[0] > 0
    assert e[-1] <= c["block_cm"], "no segment can be longer than a block"
    assert c["pools"] >= 1 and c["replicates"] >= 1 and c["block_cm"] > 0
    for cm in c["genome_cm"]:
        n = cm / c["block_cm"]
        assert np.isclose(n, round(n)) and 2 <= round(n) <= c["blocks"], \
            f"{cm} cM must be at least two whole blocks and at most the pool"
    assert c["genome_cm"] == sorted(set(c["genome_cm"]))
    assert c["diploid_samples"] >= 2 and c["haploid_ne"] > 0
    assert c["recombination_rate"] > 0 and c["mutation_rate"] > 0
    assert c["ancestry_model"] in ("smc", "smc_prime", "hudson")
    assert 0 <= c["snp_min_maf"] < 0.5 and c["snp_se_block_cm"] > 0
    assert np.isclose(c["block_cm"] / c["snp_se_block_cm"], round(c["block_cm"] / c["snp_se_block_cm"]))
    t = c["times"]
    order = ["loop_split", "loop_donor_merge", "loop_close", "ancient_admix", "left", "right", "root"]
    assert all(t[x] < t[y] for x, y in zip(order, order[1:])) and t["loop_split"] > 1
    assert c["T_max"] > t["root"]
    assert all(0 < f < 1 for f in c["fractions"].values())
    assert c["hapibd_min_output"] <= e[0]
    assert c["fit_seeds"] and c["pathfinder_paths"] > 0 and c["pathfinder_draws"] >= c["pathfinder_paths"]
    return c


def edges(c):
    return np.asarray(c["bin_edges_cm"], float)


def n_haploid(c):
    return [2 * c["diploid_samples"]] * len(POPS)


def pair_counts(c):
    return ibd.pair_counts(n_haploid(c))


# ---------------------------------------------------------------------------
# Graphs
# ---------------------------------------------------------------------------
def build(name):
    """Topology only (no parameters); the event order is the generating order."""
    d = DemographicTopology(POPS)
    if name in ("T_true", "T_alt2"):
        d.add_admixture_event("a", "aP1", "aP2")
        d.add_merge_event("aP2", "b", "bP")
        d.add_merge_event("bP", "aP1", "ab")
    else:
        d.add_merge_event("a", "b", "ab")
    if name in ("T_true", "T_alt1"):
        d.add_admixture_event("c", "cP1", "cP2")
        d.add_merge_event("ab", "cP1", "left")
        d.add_merge_event("d", "cP2", "right")
        d.add_merge_event("left", "right", "root")
    elif name == "T_alt2":
        d.add_merge_event("ab", "d", "abd")
        d.add_merge_event("abd", "c", "root")
    else:
        raise ValueError(name)
    return d


def true_topology(c):
    """The generating graph with every time, fraction and size set."""
    d = build("T_true")
    t = c["times"]; f = c["fractions"]
    d.set_uniform_ne(int(c["haploid_ne"]))
    d.set_admixture_parameters("a", t["loop_split"], f["loop"], "aP1")
    d.set_merge_time("bP", t["loop_donor_merge"])
    d.set_merge_time("ab", t["loop_close"])
    d.set_admixture_parameters("c", t["ancient_admix"], f["ancient"], "cP1")
    d.set_merge_time("left", t["left"])
    d.set_merge_time("right", t["right"])
    d.set_merge_time("root", t["root"])
    d.finalize_root()
    d.is_valid()
    return d


def msprime_demography(c):
    return to_msprime(true_topology(c))


def candidates():
    out = []
    for name in CANDIDATES:
        d = build(name)
        out.append(dict(id=name, events=d.ordered_events, nodes=list(d.nodes), n_admixture=d.n_admix))
    return out


def replay(events):
    return _replay(POPS, events)


# Which event of each candidate stands for which true parameter.  Only events
# with the same meaning in the candidate and the truth are mapped: in T_alt1 the
# a-b merge closes the (omitted) loop, so it is compared with loop_close; in
# T_alt2 the root joins c to (ab, d), which has no counterpart in the truth.
SEMANTICS = {
    "T_true": {"a": "loop_split", "bP": "loop_donor_merge", "ab": "loop_close", "c": "ancient_admix",
               "left": "left", "right": "right", "root": "root"},
    "T_alt1": {"ab": "loop_close", "c": "ancient_admix", "left": "left", "right": "right", "root": "root"},
    "T_alt2": {"a": "loop_split", "bP": "loop_donor_merge", "ab": "loop_close"},
}
FRACTIONS = {"a": ("loop_fraction", "loop"), "c": ("ancient_fraction", "ancient")}


def event_parameters(candidate, c):
    """[{name, variable, index, truth}] for every mapped time and fraction.
    Fractions are to the FIRST listed source (aP1, cP1), as in the Stan models."""
    names = SEMANTICS[candidate["id"]]
    out = []
    admixture_index = 0
    for k, ev in enumerate(candidate["events"]):
        node = ev["parent"] if ev["type"] == "MERGE" else ev["child"]
        if node in names:
            out.append(dict(name=names[node], variable="cumulative_times", index=k,
                            truth=c["times"][names[node]]))
        if ev["type"] == "ADMIXTURE":
            if node in names:
                label, key = FRACTIONS[node]
                out.append(dict(name=label, variable="admixture_fractions", index=admixture_index,
                                truth=c["fractions"][key]))
            admixture_index += 1
    return out


# ---------------------------------------------------------------------------
# Subsamples: every genome length draws its own block sets from one pool
# ---------------------------------------------------------------------------
def selections(c, pool, cm):
    return blocks.selections(c["seed"], pool, c["blocks"], round(cm / c["block_cm"]), c["replicates"], int(cm))


def tasks(c, sources, models, replicate_limit=None):
    """Every (cm, rep, source, model, candidate) fit, in a fixed order that the
    workers slice.  SNP-only never reads IBD, so it is fitted once per
    subsample with source "none" instead of once per IBD source."""
    reps = range(min(c["replicates"], replicate_limit or c["replicates"]))
    out = []
    for cm in c["genome_cm"]:
        for rep in reps:
            for model in models:
                for source in (["none"] if model == "snp" else sources):
                    for cand in CANDIDATES:
                        out.append((cm, rep, source, model, cand))
    return out


# ---------------------------------------------------------------------------
# Observations, with this study's configuration bound in
# ---------------------------------------------------------------------------
def true_ibd(ts, c):
    return ibd.true_ibd(ts, POPS, edges(c), c["recombination_rate"])


def hapibd(ts, c, folder, command):
    return ibd.hapibd(ts, POPS, edges(c), c["recombination_rate"], folder, command,
                      c["hapibd_min_seed"], c["hapibd_min_output"], c["hapibd_threads"])


def snp_summaries(ts, c):
    return snp.snp_summaries(ts, POPS, c["block_cm"], c["snp_se_block_cm"], c["recombination_rate"],
                             c["snp_min_maf"], c["times"]["root"] if c["snp_ancestral_only"] else None)


def aggregate(block_list, chosen, source, c):
    """Observed data; SNP-only fits ("none") read the true-IBD slot, which they ignore."""
    return blocks.aggregate(block_list, chosen, "true" if source == "none" else source,
                            c["block_cm"], n_haploid(c))


def stan_data(candidate, obs, c, model):
    return stan.stan_data(replay(candidate["events"]), obs, MODEL_NAMES[model],
                          edges=edges(c), n_haploid=n_haploid(c), T_max=c["T_max"])
