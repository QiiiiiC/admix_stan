"""Admixture-graph DSL.

A `DemographicTopology` is built backwards in time from its leaves by two events:

  add_merge_event(c1, c2, parent)      two lineages join into one ancestor
  add_admixture_event(child, p1, p2)   one lineage splits into two sources

The order events are added IS their temporal order: the Stan models integrate
`times` as positive increments (cumulative_times = cumulative_sum(times)), so
event k of `ordered_events` is cumulative_times[k].

Conventions every model and study relies on:
  * Ne on a node is DIPLOID (`set_node_ne(name, haploid // 2)`); the Stan models
    use HAPLOID size, so true haploid Ne = 2 * node.ne.
  * Stan node index a (1-based) == position of the node in `nodes` insertion
    order, so `Ne[:, k]` <-> `list(dem.nodes)[k]`.
  * MERGE event time = nodes[parent].time_start; ADMIXTURE event time =
    nodes[child].time_end.
  * Stan admixture_fractions[i] is the fraction from the FIRST listed parent of
    the i-th (chronological) admixture event.

Plotting lives with the studies, not here.
"""
from __future__ import annotations

import itertools

import numpy as np


class Node:
    """A population branch.  Ancestors are further in the past, descendants closer to the present."""

    def __init__(self, name):
        self.name = name
        self.descendants = []
        self.ancestors = []
        self.admixture_fractions = {}
        self.ne = None          # DIPLOID effective size
        self.time_start = None  # generations ago this branch appears (young end)
        self.time_end = None    # generations ago it merges/splits (old end)

    def add_descendant(self, node):
        self.descendants.append(node)

    def add_ancestor(self, node, fraction=1.0):
        self.ancestors.append(node)
        self.admixture_fractions[node.name] = fraction

    @property
    def is_root(self):
        return not self.ancestors

    @property
    def is_leaf(self):
        return not self.descendants

    @property
    def is_admixed(self):
        return len(self.ancestors) > 1


class DemographicTopology:
    def __init__(self, leaf_names):
        self.nodes = {}
        self.ordered_events = []
        self.initial_leaves = list(leaf_names)
        self.n_admix = 0
        for name in leaf_names:
            node = Node(name)
            node.time_start = 0.0
            self.nodes[name] = node

    def get_node(self, name):
        if name not in self.nodes:
            raise ValueError(f"Node '{name}' not found.")
        return self.nodes[name]

    # ---- topology -----------------------------------------------------------
    def add_merge_event(self, child_1_name, child_2_name, new_parent_name):
        if new_parent_name in self.nodes:
            raise ValueError(f"Node '{new_parent_name}' already exists.")
        children = [self.get_node(child_1_name), self.get_node(child_2_name)]
        parent = Node(new_parent_name)
        for child in children:
            child.add_ancestor(parent, fraction=1.0)
            parent.add_descendant(child)
        self.nodes[new_parent_name] = parent
        self.ordered_events.append({"id": len(self.ordered_events) + 1, "type": "MERGE",
                                    "children": [child_1_name, child_2_name],
                                    "parent": new_parent_name})

    def add_admixture_event(self, child_name, parent_1_name, parent_2_name):
        if parent_1_name in self.nodes or parent_2_name in self.nodes:
            raise ValueError("New parent nodes must have unique names.")
        child = self.get_node(child_name)
        for name in (parent_1_name, parent_2_name):
            parent = Node(name)
            child.add_ancestor(parent, fraction=None)   # set by set_admixture_parameters
            parent.add_descendant(child)
            self.nodes[name] = parent
        self.ordered_events.append({"id": len(self.ordered_events) + 1, "type": "ADMIXTURE",
                                    "child": child_name,
                                    "parents": [parent_1_name, parent_2_name]})
        self.n_admix += 1

    # ---- parameters ---------------------------------------------------------
    def set_node_ne(self, node_name, ne):
        """DIPLOID effective size."""
        self.get_node(node_name).ne = ne

    def set_uniform_ne(self, ne_haploid):
        for name in self.nodes:
            self.set_node_ne(name, ne_haploid // 2)

    def set_merge_time(self, parent_name, time):
        parent = self.get_node(parent_name)
        parent.time_start = time
        for child in parent.descendants:
            child.time_end = time

    def set_admixture_parameters(self, child_name, time, fraction_parent_1, parent_1_name):
        child = self.get_node(child_name)
        child.time_end = time
        if parent_1_name not in [a.name for a in child.ancestors]:
            raise ValueError(f"Parent {parent_1_name} is not an ancestor of {child_name}")
        for ancestor in child.ancestors:
            ancestor.time_start = time
            child.admixture_fractions[ancestor.name] = (
                fraction_parent_1 if ancestor.name == parent_1_name else 1.0 - fraction_parent_1)

    def finalize_root(self):
        for node in self.nodes.values():
            if node.is_root:
                node.time_end = float("inf")

    # ---- checks and matrix form ---------------------------------------------
    def is_valid(self):
        """Events consume active lineages, end at one root, and times are ordered."""
        alive = set(self.initial_leaves)
        for ev in self.ordered_events:
            consumed = ev["children"] if ev["type"] == "MERGE" else [ev["child"]]
            if not set(consumed) <= alive:
                raise ValueError(f"Event {ev['id']}: {consumed} are not active lineages.")
            alive -= set(consumed)
            alive |= {ev["parent"]} if ev["type"] == "MERGE" else set(ev["parents"])
        if len(alive) != 1:
            raise ValueError(f"Ends with {len(alive)} roots ({alive}); must be exactly 1.")
        for node in self.nodes.values():
            if (node.time_start is not None and node.time_end is not None
                    and node.time_start > node.time_end):
                raise ValueError(f"{node.name} starts at {node.time_start} after it ends at {node.time_end}.")
        return True

    def get_topology_matrix_representation(self):
        """(migration_matrices, event_types, admixture_map, admixture_map_id).

        One n_nodes x n_nodes matrix per event.  MERGE: children map to the parent,
        everything else to itself.  ADMIXTURE: identity -- the fractions are
        parameters, so the Stan model builds that matrix itself from admixture_map.
        """
        index = {name: i for i, name in enumerate(self.nodes)}
        n = len(index)
        matrices, types, admix, admix_id = [], [], {}, {}
        for k, ev in enumerate(self.ordered_events):
            types.append(ev["type"])
            mat = np.eye(n)
            if ev["type"] == "MERGE":
                for child in ev["children"]:
                    mat[index[child], index[child]] = 0
                    mat[index[child], index[ev["parent"]]] = 1
            else:
                admix[k] = {"child": ev["child"], "parents": ev["parents"]}
                admix_id[k] = {"child": index[ev["child"]],
                               "parents": [index[p] for p in ev["parents"]]}
            matrices.append(mat)
        return matrices, types, admix, admix_id

    def event_times(self):
        """True time of every event, in event order (see module conventions)."""
        return [self.nodes[ev["parent"]].time_start if ev["type"] == "MERGE"
                else self.nodes[ev["child"]].time_end for ev in self.ordered_events]

    def admixture_fractions(self):
        """True fraction from the FIRST listed parent, per admixture event in order."""
        return [self.nodes[ev["child"]].admixture_fractions[ev["parents"][0]]
                for ev in self.ordered_events if ev["type"] == "ADMIXTURE"]

    def summary(self):
        rows = [f"{'node':<12} {'start':>8} {'end':>8} {'Ne(dip)':>9}  ancestors"]
        for n in sorted(self.nodes.values(), key=lambda x: -1 if x.time_start is None else x.time_start):
            anc = ", ".join(f"{a.name}({n.admixture_fractions.get(a.name)})" for a in n.ancestors) or "ROOT"
            rows.append(f"{n.name:<12} {str(n.time_start):>8} {str(n.time_end):>8} {str(n.ne):>9}  {anc}")
        return "\n".join(rows)


# ---------------------------------------------------------------------------
# Working with event lists
# ---------------------------------------------------------------------------
def replay(leaves, events):
    """Rebuild a topology from a saved `ordered_events` list."""
    d = DemographicTopology(leaves)
    for ev in events:
        if ev["type"] == "MERGE":
            d.add_merge_event(*ev["children"], ev["parent"])
        else:
            d.add_admixture_event(ev["child"], *ev["parents"])
    return d


def valid_event_orders(leaves, events):
    """Every permutation of `events` that is a valid backwards-in-time history.

    The same graph can be written with different event orders (e.g. which of two
    independent merges is older); each order is a different Stan model because
    times are increments.  An order is valid if every event consumes lineages
    that are alive at that point and it ends at a single root.
    """
    valid = []
    for order in itertools.permutations(events):
        alive = set(leaves)
        for ev in order:
            consumed = set(ev["children"] if ev["type"] == "MERGE" else [ev["child"]])
            if not consumed <= alive:
                break
            alive -= consumed
            alive.update([ev["parent"]] if ev["type"] == "MERGE" else ev["parents"])
        else:
            if len(alive) == 1:
                valid.append(list(order))
    return valid


def clades(leaves, events):
    """Leaf set under every node.  Admixture sources start a fresh clade named
    after themselves, so a node's clade says which lineages it gathers."""
    out = {name: {name} for name in leaves}
    for ev in events:
        if ev["type"] == "ADMIXTURE":
            for p in ev["parents"]:
                out[p] = {p}
        else:
            out[ev["parent"]] = set.union(*(out[c] for c in ev["children"]))
    return out
