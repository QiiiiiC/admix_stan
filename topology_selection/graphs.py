"""Admixture graphs as topologies, their canonical keys, and the moves between them.

A graph is a set of events.  Each event consumes branches and creates branches:

    merge  (x, y) -> (p,)          admix  (x,) -> (y, z)      y, z: the two sources

Leaves are branches with no creator; the root is the one branch with no consumer.
There are no heights and no event order: this is exactly what the per-branch model
(branch_model/) takes, and `to_topology` hands a graph to it.

Moves (README.md, section 2):
    re-merge   detach a branch from its merge, merge it with another branch
    add        split a branch by a new admixture; its second source merges with a branch
    remove     detach one source of an admixture and fuse the other back
    relocate   remove then add, as one move (only the end point is scored): moves an
               admixture to another population without passing through a graph with
               one admixture fewer
"""
from __future__ import annotations

from dataclasses import dataclass, field
import itertools
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "branch_model")]
from methods.demography import DemographicTopology


@dataclass
class Graph:
    leaves: tuple
    events: dict = field(default_factory=dict)      # event id -> (kind, inputs, outputs)
    counter: int = 0                                # for fresh branch and event names

    def copy(self):
        return Graph(self.leaves, dict(self.events), self.counter)

    def fresh(self, prefix):
        self.counter += 1
        return f"{prefix}{self.counter}"

    # ---- structure -------------------------------------------------------
    def creator(self):
        return {b: e for e, (_, _, outs) in self.events.items() for b in outs}

    def consumer(self):
        return {b: e for e, (_, ins, _) in self.events.items() for b in ins}

    def branches(self):
        return list(self.leaves) + [b for _, _, outs in self.events.values() for b in outs]

    def root(self):
        con = self.consumer()
        (r,) = [b for b in self.branches() if b not in con]
        return r

    def down(self, b):
        """b and every branch that flows into it (its subnetwork towards the present)."""
        cre, out, todo = self.creator(), set(), [b]
        while todo:
            x = todo.pop()
            if x not in out:
                out.add(x)
                if x in cre:
                    todo.extend(self.events[cre[x]][1])
        return out

    def n_admix(self):
        return sum(kind == "admix" for kind, _, _ in self.events.values())

    def listing(self):
        """The events in an order that respects "must come after" (topological sort)."""
        cre, done, order = self.creator(), set(), []
        while len(order) < len(self.events):
            ready = sorted(e for e in self.events if e not in done and
                           all(b not in cre or cre[b] in done for b in self.events[e][1]))
            if not ready:
                raise ValueError("the graph has a cycle")
            done.update(ready); order += ready
        return order

    # ---- primitive operations --------------------------------------------
    def _replace_input(self, old, new):
        e = self.consumer().get(old)
        if e is not None:
            kind, ins, outs = self.events[e]
            self.events[e] = (kind, tuple(new if b == old else b for b in ins), outs)

    def detach(self, b):
        """b is consumed by merge(b, c) -> p: remove that merge; c takes p's place.
        Returns c."""
        e = self.consumer()[b]
        kind, ins, (p,) = self.events[e]
        assert kind == "merge"
        c = ins[1] if ins[0] == b else ins[0]
        del self.events[e]
        self._replace_input(p, c)
        return c

    def attach(self, b, t):
        """New merge(b, t) -> t'; t' takes over t's consumer (t' is the new root if t was)."""
        t2 = self.fresh("m")
        self._replace_input(t, t2)
        self.events[self.fresh("E")] = ("merge", (b, t), (t2,))

    def split(self, x):
        """New admixture x -> (y, z); y takes over x's consumer, z is left dangling."""
        y, z = self.fresh("s"), self.fresh("s")
        self._replace_input(x, y)
        self.events[self.fresh("A")] = ("admix", (x,), (y, z))
        return y, z

    def unsplit(self, a, keep):
        """Remove admixture a = x -> (keep, dropped) whose dropped source is already detached."""
        _, (x,), _ = self.events.pop(a)
        self._replace_input(keep, x)

    # ---- identity ----------------------------------------------------------
    def key(self):
        """Canonical string.  A branch is labelled by its subnetwork: a leaf by its name,
        a merge output by its sorted inputs, an admixture source by [input] and a side.
        The two sources of an admixture are exchangeable (f <-> 1 - f), so the key is the
        minimum over the 2^k side assignments."""
        admixes = sorted(e for e, (kind, _, _) in self.events.items() if kind == "admix")
        cre, root, best = self.creator(), self.root(), None
        for flips in itertools.product((0, 1), repeat=len(admixes)):
            flip, memo = dict(zip(admixes, flips)), {}

            def label(b):
                if b not in memo:
                    if b not in cre:
                        memo[b] = b
                    else:
                        e = cre[b]
                        kind, ins, outs = self.events[e]
                        memo[b] = ("(" + ",".join(sorted(label(x) for x in ins)) + ")" if kind == "merge"
                                   else f"[{label(ins[0])}]{(outs.index(b) + flip[e]) % 2 + 1}")
                return memo[b]
            s = label(root)
            best = s if best is None or s < best else best
        return best


# ---------------------------------------------------------------------------
# Conversion to and from the repository's graph language
# ---------------------------------------------------------------------------
def to_topology(g):
    """A DemographicTopology whose events follow a valid listing; branch names kept."""
    d = DemographicTopology(list(g.leaves))
    for e in g.listing():
        kind, ins, outs = g.events[e]
        if kind == "merge":
            d.add_merge_event(ins[0], ins[1], outs[0])
        else:
            d.add_admixture_event(ins[0], outs[0], outs[1])
    return d


def from_topology(dem):
    g = Graph(tuple(dem.initial_leaves))
    for ev in dem.ordered_events:
        if ev["type"] == "MERGE":
            g.events[g.fresh("E")] = ("merge", tuple(ev["children"]), (ev["parent"],))
        else:
            g.events[g.fresh("A")] = ("admix", (ev["child"],), tuple(ev["parents"]))
    return g


def caterpillar(leaves):
    g = Graph(tuple(leaves))
    cur = leaves[0]
    for x in leaves[1:]:
        p = g.fresh("m")
        g.events[g.fresh("E")] = ("merge", (cur, x), (p,))
        cur = p
    return g


# ---------------------------------------------------------------------------
# The space and the moves
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Space:
    max_admix: int = 1
    bubbles: bool = False         # an admixture whose two sources merge with each other next
    ghost: bool = False           # a source that joins only after every other lineage has merged
    targets: str = "any"          # "any" or "leaves": which branches may be admixed


def allowed(g, sp):
    """Constraints a graph must satisfy.  Checked on every proposed graph; the same test
    in both directions keeps every move reversible."""
    if g.n_admix() > sp.max_admix:
        return False
    con, root = g.consumer(), g.root()
    for kind, ins, outs in g.events.values():
        if kind != "admix":
            continue
        if sp.targets == "leaves" and ins[0] not in g.leaves:
            return False
        y, z = outs
        together = con.get(y) is not None and con.get(y) == con.get(z)
        if together and g.events[con[y]][2][0] == root:
            return False              # the root lineage admixed with itself: degenerate
        if together and not sp.bubbles:
            return False
        if not sp.ghost and any(s in con and g.events[con[s]][0] == "merge" and g.events[con[s]][2][0] == root
                                for s in outs):
            return False
    return True


def moves(g, sp):
    """Every move from g as (kind, args, new graph), before the space's constraints."""
    out, root, con = [], g.root(), g.consumer()
    # re-merge: covers NNI, SPR and re-rooting of the tree part, and moving an admixture source
    for b in g.branches():
        if b not in con or g.events[con[b]][0] != "merge":
            continue
        h = g.copy()
        c = h.detach(b)
        below = h.down(b)
        for t in h.branches():
            if t not in below and t != c:
                n = h.copy(); n.attach(b, t)
                out.append(("re-merge", (b, t), n))
    # add an admixture on x; its second source z merges with t
    if g.n_admix() < sp.max_admix:
        for x in g.branches():
            if x == root or (sp.targets == "leaves" and x not in g.leaves):
                continue
            h = g.copy()
            _, z = h.split(x)
            below = h.down(z)
            for t in h.branches():
                if t not in below:
                    n = h.copy(); n.attach(z, t)
                    out.append(("add", (x, t), n))
    # remove an admixture: drop one source (it must be consumed by a merge), fuse the other back
    for a, (kind, ins, outs) in g.events.items():
        if kind != "admix":
            continue
        for drop in outs:
            if drop in con and g.events[con[drop]][0] == "merge":
                keep = outs[1] if drop == outs[0] else outs[0]
                h = g.copy(); h.detach(drop); h.unsplit(a, keep)
                out.append(("remove", (a, drop), h))
    return out


def neighbours(g, sp, relocate=False):
    """key -> [multiplicity, graph]: the graphs one valid move away.  The proposal picks
    a move uniformly, so q(g -> g') = multiplicity / total multiplicity.  `relocate`
    adds every remove-then-add path as one move (it is its own reverse)."""
    here, out = g.key(), {}
    found = [n for _, _, n in moves(g, sp)]
    if relocate:
        found += [n2 for kind, _, n in moves(g, sp) if kind == "remove"
                  for kind2, _, n2 in moves(n, add_space(sp, n)) if kind2 == "add"]
    for n in found:
        if not allowed(n, sp):
            continue
        k = n.key()
        if k != here:
            out.setdefault(k, [0, n])[0] += 1
    return out


def add_space(sp, g):
    """The add moves from an intermediate graph are not limited by its own admixture count."""
    return Space(max_admix=g.n_admix() + 1, bubbles=sp.bubbles, ghost=sp.ghost, targets=sp.targets)


def enumerate_space(sp, leaves, relocate=False):
    """Every graph reachable from one tree: key -> graph, and key -> neighbours."""
    start = caterpillar(list(leaves))
    seen, nbrs, frontier = {start.key(): start}, {}, [start]
    while frontier:
        nxt = []
        for g in frontier:
            nbrs[g.key()] = neighbours(g, sp, relocate)
            for k, (_, h) in nbrs[g.key()].items():
                if k not in seen:
                    seen[k] = h; nxt.append(h)
        frontier = nxt
    return seen, nbrs
