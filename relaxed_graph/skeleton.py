"""The skeleton [admixture, merge, ..., merge, root] on L leaves (MATHS.md section 2).

Demes 1..L are the leaves (0-based here: 0..L-1), deme L the ghost (the admixture's
second source).  A vertex v = (x, p_1, ..., p_{L-1}) picks the admixed leaf x and,
for each merge slot, a pair p = (a, b), a < b, whose deme b moves into deme a.  The
last slot is the compulsory root.  Vertices are discrete graphs; several vertices
can be the same graph: the admixture's two sources are exchangeable (swap them and
f <-> 1 - f), and a merge that involves an empty deme does nothing, so the root then
joins three or more lineages at once.
"""
from __future__ import annotations

from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from methods.demography import DemographicTopology


def _canon(c):
    if isinstance(c, tuple):
        return f"{c[0]}{c[1]}"
    if isinstance(c, frozenset):
        return "(" + ",".join(sorted(_canon(y) for y in c)) + ")"
    return str(c)


class Skeleton:
    def __init__(self, leaves):
        self.leaves = list(leaves)
        self.L = len(self.leaves)
        self.D = self.L + 1
        self.n_merge = self.L - 1                       # L+1 lineages -> 2 before the root
        self.K = self.n_merge + 2
        self.pairs = list(combinations(range(self.D), 2))
        self.vertices = [(x, *ps) for x in range(self.L)
                         for ps in product(range(len(self.pairs)), repeat=self.n_merge)]

    # ------------------------------------------------------------------
    def replay(self, v, swap=False):
        """Walk the slots.  Returns (events, contents): events lists the effective
        events as (kind, descriptors); contents the final deme descriptors.  `swap`
        names the staying source 2 and the ghost 1."""
        x, *ps = v
        content = {d: self.leaves[d] for d in range(self.L)}
        content[self.L] = None                          # ghost empty before the admixture
        leaf = self.leaves[x]
        events = [("admix", leaf)]
        content[x], content[self.L] = ((leaf, 2), (leaf, 1)) if swap else ((leaf, 1), (leaf, 2))
        for p in ps:
            a, b = self.pairs[p]
            if content[a] is None or content[b] is None:
                continue                                # no lineage moves: a no-op slot
            events.append(("merge", frozenset([content[a], content[b]])))
            content[a], content[b] = frozenset([content[a], content[b]]), None
        live = [c for c in content.values() if c is not None]
        events.append(("root", frozenset(live)))
        return events, live

    def key(self, v):
        """Graph identity: the ordered effective events in terms of lineage content,
        independent of deme labels and of which admixture source is called 1.  A
        string with sorted sub-terms, so equal graphs give equal keys in any process
        (frozenset iteration order is not canonical)."""
        return min(" ; ".join(f"{kind} {_canon(what)}" for kind, what in self.replay(v, swap)[0])
                   for swap in (False, True))

    def is_binary(self, v):
        return len(self.replay(v)[1]) == 2

    def label(self, v):
        def name(c):
            if isinstance(c, tuple):
                return f"{c[0]}{c[1]}"
            if isinstance(c, frozenset):
                return "(" + ",".join(sorted(name(y) for y in c)) + ")"
            return str(c)
        out = []
        for kind, what in self.replay(v)[0]:
            if kind == "admix":
                out.append(f"admix {what}")
            else:
                out.append(("root " if kind == "root" else "") + "+".join(sorted(name(y) for y in what)))
        return "; ".join(out)

    def graphs(self):
        """Distinct graphs: key -> list of vertex indices (the prior mass of a graph is
        its share of vertices)."""
        out = {}
        for i, v in enumerate(self.vertices):
            out.setdefault(self.key(v), []).append(i)
        return out

    # ------------------------------------------------------------------
    def matrices(self, w_admix, w_merge, f):
        """Numpy twin of slot_matrices() in stan/graph.stanfunctions (eq. 9-11)."""
        D, L = self.D, self.L
        M = [np.eye(D)]
        for x in range(L):
            M[0][x, x] = 1 - w_admix[x] * (1 - f)
            M[0][x, L] = w_admix[x] * (1 - f)
        for k in range(self.n_merge):
            m = np.eye(D)
            for p, (a, b) in enumerate(self.pairs):
                m[b, b] -= w_merge[k][p]
                m[b, a] += w_merge[k][p]
            M.append(m)
        root = np.zeros((D, D)); root[:, 0] = 1
        M.append(root)
        return M

    def one_hot(self, v):
        x, *ps = v
        return np.eye(self.L)[x], [np.eye(len(self.pairs))[p] for p in ps]

    def topology(self, v):
        """The same graph in the repository's DSL (binary vertices only), with node
        names that follow the deme contents; used to test against methods' models."""
        if not self.is_binary(v):
            raise ValueError("a vertex whose root joins more than two lineages has no binary DSL graph")
        x, *ps = v
        d = DemographicTopology(self.leaves)
        node = {k: self.leaves[k] for k in range(self.L)}
        leaf = self.leaves[x]
        node[x], node[self.L] = f"{leaf}_1", f"{leaf}_2"
        d.add_admixture_event(leaf, node[x], node[self.L])
        live = set(range(self.D))
        for k, p in enumerate(ps):
            a, b = self.pairs[p]
            if a not in live or b not in live:
                continue
            parent = f"m{k + 1}"
            d.add_merge_event(node[a], node[b], parent)
            node[a] = parent; live.discard(b)
        a, b = sorted(live)
        d.add_merge_event(node[a], node[b], "root")
        return d
