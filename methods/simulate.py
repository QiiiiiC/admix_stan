"""msprime simulation of a generating history.

Sizes are given HAPLOID everywhere in this module and halved on the way into
msprime (which takes diploid sizes).  A branch may grow exponentially between
its young end t0 and its old end t1; see `add_branch`.
"""
from __future__ import annotations

import numpy as np


def growth_rate(young, old, t0, t1):
    """msprime rate (per generation, backwards) taking size `young` at t0 to `old` at t1."""
    if young == old or not np.isfinite(t1):
        return 0.0
    return float(np.log(young / old) / (t1 - t0))


def size_at(young, old, t0, t1, t):
    """Haploid size at time t of a branch exponential between (t0, young) and (t1, old)."""
    t = np.asarray(t, float)
    if young == old or not np.isfinite(t1):
        return np.full(t.shape, float(young))
    return young * np.exp(np.log(old / young) * np.clip((t - t0) / (t1 - t0), 0, 1))


def add_branch(demography, name, young, old, t0, t1):
    """Add one population living on [t0, t1] with haploid size young -> old.

    msprime measures growth from time 0, so an ancestral branch (t0 > 0) gets a
    population_parameters_change AT t0 setting its young-end size and rate; that
    anchors size(t0) exactly.  Growth is switched off at t1 so no branch keeps
    changing after its lineages have moved on.
    """
    rate = growth_rate(young, old, t0, t1)
    demography.add_population(name=name, initial_size=young / 2,
                              growth_rate=rate if t0 == 0 else 0.0)
    if rate and t0 > 0:
        demography.add_population_parameters_change(time=t0, population=name,
                                                    initial_size=young / 2, growth_rate=rate)
    if rate:
        demography.add_population_parameters_change(time=t1, population=name, growth_rate=0.0)


def event_time(dem, ev):
    return dem.nodes[ev["parent"]].time_start if ev["type"] == "MERGE" else dem.nodes[ev["child"]].time_end


def to_msprime(dem, sizes=None):
    """msprime.Demography for a fully parameterised DemographicTopology.

    sizes: optional {node: (young, old)} HAPLOID sizes for branches that change
    size; every other node is constant at 2 * node.ne.  Node names must be valid
    msprime identifiers (no dots) for a generating history.
    """
    import msprime
    sizes = sizes or {}
    d = msprime.Demography()
    for name, node in dem.nodes.items():
        if name in sizes:
            young, old = sizes[name]
        else:
            if node.ne is None:
                raise ValueError(f"Node {name} has no Ne set")
            young = old = 2.0 * node.ne
        t1 = node.time_end if node.time_end is not None else np.inf
        add_branch(d, name, young, old, node.time_start, t1)
    for ev in dem.ordered_events:
        t = event_time(dem, ev)
        if t is None:
            raise ValueError(f"Time not set for event {ev}")
        if ev["type"] == "MERGE":
            d.add_population_split(time=t, derived=list(ev["children"]), ancestral=ev["parent"])
        else:
            fractions = dem.nodes[ev["child"]].admixture_fractions
            p1, p2 = ev["parents"]
            if fractions.get(p1) is None:
                raise ValueError(f"Admixture fractions not set for {ev['child']}")
            d.add_admixture(time=t, derived=ev["child"], ancestral=[p1, p2],
                            proportions=[fractions[p1], fractions[p2]])
    d.sort_events()
    d.validate()
    return d


def simulate_block(demography, samples, block_cm, recombination_rate, mutation_rate, seed,
                   model="smc"):
    """One independent genome block: ancestry + mutations on a discrete genome.

    samples: {population: number of DIPLOID individuals}.  Block length in bp is
    block_cm / (100 * recombination_rate), i.e. a uniform map.
    """
    import msprime
    ancestry = msprime.sim_ancestry(
        samples=samples, ploidy=2, demography=demography,
        sequence_length=round(block_cm / (100 * recombination_rate)),
        recombination_rate=recombination_rate, random_seed=seed,
        model=model, discrete_genome=True)
    return msprime.sim_mutations(ancestry, rate=mutation_rate, random_seed=seed, discrete_genome=True)


def block_seed(*key):
    """Deterministic, never-zero 32-bit seed from an integer key (seed, scenario, pool, block, ...)."""
    return int(np.random.SeedSequence(list(key)).generate_state(1)[0]) or 1
