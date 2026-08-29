"""Every rooted admixture graph on three leaves with at most one admixture event.

RESULT: 3 trees + 18 one-admixture graphs = 21 distinct models.

------------------------------------------------------------------------------
Counting, and where the natural enumeration goes wrong
------------------------------------------------------------------------------
No admixture: 3 rooted trees on 3 labelled leaves -- ((A,B),C), ((A,C),B),
((B,C),A).  Nothing subtle.

One admixture, admixed leaf X: going backwards X splits into two source
lineages X.1 and X.2, so the graph is a rooted binary tree on the FOUR lineages
{X.1, X.2, A, B}.  There are (2n-3)!! = 15 such trees.  Two corrections:

  (i)  drop the trees where X.1 and X.2 are siblings.  There they coalesce with
       each other before meeting anything else, so no lineage ever separates
       them and the admixture fraction is unidentifiable -- it is just the tree
       ((A,B),X) with a decorative extra node.  That removes 3.

  (ii) quotient by the swap X.1 <-> X.2.  Relabelling the two sources and
       sending f -> 1-f is the SAME graph.  The swap has no fixed points among
       the surviving 12 (a fixed point would need X.1 and X.2 in symmetric
       positions, i.e. siblings, already removed), so the 12 pair up exactly.

  12 / 2 = 6 per admixed leaf, and 3 x 6 = 18.

The obvious hand enumeration -- "the first merge joins an admixture branch to
one of the other two leaves (2 choices), leaving three lineages (3 topologies)",
giving 3 x 2 x 3 = 18 -- lands on the right TOTAL by cancelling two mistakes:

  * it DOUBLE-COUNTS the balanced graph.  ((X.1,A),(X.2,B)) and ((X.1,B),(X.2,A))
    appear as separate entries, but they are the same graph read with the source
    labels swapped.  -3 (one per admixed leaf).

  * it MISSES the graph where the two non-admixed leaves merge FIRST:
    (((A,B),X.1),X.2).  Nothing in the first merge involves an admixture branch,
    so the stated rule forbids it -- yet it is the most common real scenario
    there is: a population that is part local clade, part deep outgroup.  +3.

15 + 3 = 18 either way, but the SET differs by three graphs.  All 18 are built
here; `outside_first_merge_rule` marks the three the hand rule would have missed.

------------------------------------------------------------------------------
Event ordering
------------------------------------------------------------------------------
The Stan model integrates `times` as positive INCREMENTS (cumulative_times =
cumulative_sum(times)), so the order events are added IS their temporal order.
Two conventions are imposed:

  * the admixture event is always event 1.  It sits on a leaf branch, so it is
    the most recent event on that lineage; the only thing this forbids is an
    admixture older than a merge between the other two leaves.
  * merges follow in post-order, so every merge comes after both of its
    children's events.  This is forced, not a convention.

Leaf ORDER (`DemographicTopology(leaves)`) is always the saved pop_order, since
that is what indexes the data matrix rows.  It carries no topological meaning --
the events do.
"""
import os
import sys
from itertools import combinations

_RD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # real_data/
_NP = os.path.dirname(_RD)                                          # new_pipeline/
for _d in (_RD, _NP):
    if _d not in sys.path:
        sys.path.insert(0, _d)
from demography import DemographicTopology  # noqa: E402


# --------------------------------------------------------------------------
# rooted binary trees on a labelled set, as nested 2-tuples
# --------------------------------------------------------------------------
def rooted_trees(leaves):
    """All rooted binary trees on `leaves`, as nested tuples.  (2n-3)!! of them.

    Each bipartition is generated once by anchoring the first element in the
    left part, which is what stops (A,B) and (B,A) both appearing.
    """
    leaves = list(leaves)
    if len(leaves) == 1:
        return [leaves[0]]
    head, rest = leaves[0], leaves[1:]
    out = []
    for k in range(len(rest) + 1):
        for combo in combinations(rest, k):
            left = [head] + list(combo)
            right = [x for x in rest if x not in combo]
            if not right:
                continue
            for lt in rooted_trees(left):
                for rt in rooted_trees(right):
                    out.append((lt, rt))
    return out


def canonical(tree):
    """Order-independent string form, so two builds of one shape compare equal."""
    if isinstance(tree, str):
        return tree
    a, b = canonical(tree[0]), canonical(tree[1])
    return "(" + ",".join(sorted([a, b])) + ")"


def relabel(tree, mapping):
    if isinstance(tree, str):
        return mapping.get(tree, tree)
    return (relabel(tree[0], mapping), relabel(tree[1], mapping))


def are_siblings(tree, x, y):
    if isinstance(tree, str):
        return False
    a, b = tree
    if {a, b} == {x, y}:
        return True
    return are_siblings(a, x, y) or are_siblings(b, x, y)


def first_merge_involves(tree, sources):
    """Does the earliest merge (a cherry deepest in the future) touch a source?

    Only used to flag the three graphs the hand rule would drop.  A cherry of
    two plain leaves that is not the whole tree means two non-admixed leaves
    merged first.
    """
    def cherries(t):
        if isinstance(t, str):
            return []
        a, b = t
        here = [{a, b}] if isinstance(a, str) and isinstance(b, str) else []
        return here + cherries(a) + cherries(b)
    return all(c & set(sources) for c in cherries(tree))


# --------------------------------------------------------------------------
# tree -> DemographicTopology
# --------------------------------------------------------------------------
def _emit_merges(tree, dem, counter, src_of):
    """Post-order: a merge is added only once both its children exist."""
    if isinstance(tree, str):
        return tree
    left = _emit_merges(tree[0], dem, counter, src_of)
    right = _emit_merges(tree[1], dem, counter, src_of)
    counter[0] += 1
    name = f"n{counter[0]}"
    dem.add_merge_event(left, right, name)
    return name


def build(pop_order, tree, admixed=None, sources=()):
    """Assemble the DemographicTopology.  Renames the final node to 'root'."""
    dem = DemographicTopology(list(pop_order))
    if admixed is not None:
        dem.add_admixture_event(admixed, sources[0], sources[1])
    counter = [0]
    last = _emit_merges(tree, dem, counter, sources)
    # The last merge is the root; rename it so downstream reporting can find it.
    node = dem.nodes.pop(last)
    node.name = "root"
    dem.nodes["root"] = node
    for ev in dem.ordered_events:
        if ev.get("parent") == last:
            ev["parent"] = "root"
        if ev.get("children"):
            ev["children"] = [("root" if c == last else c) for c in ev["children"]]
    return dem


def enumerate_all(pop_order):
    """All 21 models, as dicts ready for fitting.

    Returns a list of {name, newick, n_admix, admixed, builder, ...}, ordered
    trees-first then by admixed leaf, so folder indices are stable across runs.
    """
    pops = list(pop_order)
    models = []

    # ---- 3 trees -------------------------------------------------------
    for t in rooted_trees(pops):
        models.append({
            "n_admix": 0,
            "admixed": None,
            "tree": t,
            "newick": canonical(t),
            "sources": (),
            "outside_first_merge_rule": False,
        })
    models.sort(key=lambda m: m["newick"])

    # ---- 18 one-admixture graphs ---------------------------------------
    admix_models = []
    for adm in pops:
        s1, s2 = f"{adm}.1", f"{adm}.2"
        others = [p for p in pops if p != adm]
        seen = {}
        for t in rooted_trees([s1, s2] + others):
            if are_siblings(t, s1, s2):
                continue                              # (i) unidentifiable f
            swapped = relabel(t, {s1: s2, s2: s1})    # (ii) quotient by the swap
            key = min(canonical(t), canonical(swapped))
            if key in seen:
                continue
            seen[key] = True
            admix_models.append({
                "n_admix": 1,
                "admixed": adm,
                "tree": t,
                "newick": key,
                "sources": (s1, s2),
                "outside_first_merge_rule": not first_merge_involves(t, (s1, s2)),
            })
    admix_models.sort(key=lambda m: (pops.index(m["admixed"]), m["newick"]))
    models += admix_models

    for i, m in enumerate(models, 1):
        m["index"] = i
        m["name"] = f"{i:02d}_{m['newick']}"
    return models


def make_dem(pop_order, m):
    return build(pop_order, m["tree"], m["admixed"], m["sources"])


if __name__ == "__main__":
    POPS = sys.argv[1:4] if len(sys.argv) > 3 else ["GBR", "IBS", "TSI"]
    ms = enumerate_all(POPS)
    n_tree = sum(1 for m in ms if m["n_admix"] == 0)
    print(f"{len(ms)} models on {POPS}: {n_tree} trees + {len(ms)-n_tree} one-admixture\n")
    print(f"{'#':>3} {'admixed':>8} {'newick':<34} {'events':>7}  note")
    import io, contextlib
    for m in ms:
        with contextlib.redirect_stdout(io.StringIO()):   # mute demography's prints
            dem = make_dem(POPS, m)
        note = "MISSED by the 'first merge must be admixed' rule" \
            if m["outside_first_merge_rule"] else ""
        print(f"{m['index']:>3} {str(m['admixed']):>8} {m['newick']:<34} "
              f"{len(dem.ordered_events):>7}  {note}")
