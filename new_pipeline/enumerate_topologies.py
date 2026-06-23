"""
Enumerate all distinct admixture-graph topologies on a fixed set of leaf
populations, with at most `max_admix` admixture events, for the IBD / mixed
demographic model.

Design (see discussion):
  * Generate via the BACKWARD coalescent process (start with the leaves as
    active lineages; repeatedly either MERGE two active lineages into one
    ancestor, or ADMIX one active lineage into two fresh parent lineages).
    Stop when a single active lineage (the root) remains.
  * Two event-histories that yield the same DAG are the SAME topology -- the
    factorial blow-up the user noticed is just event-ordering. We therefore
    deduplicate by graph isomorphism (leaves labelled, internal nodes
    unlabelled, MERGE children unordered, ADMIX parents unordered).
  * We do NOT contract degree-2 "passthrough" admixture-source segments:
    in the IBD/demographic model every branch segment carries its own
    Ne / time, so such segments are real model structure, not f-statistic
    degeneracy.

Pruning gate (structural + counting only, as chosen):
  * structural: drop the trivial "bubble" admixture (the two admixture
    sources merge directly back with each other and nothing else) -- the
    admixture proportion is non-identifiable there.
  * counting: keep only graphs with #params p <= #independent statistics m.
    For the IBD model m = (n*(n+1)/2) * n_bins, which is large, so this
    rarely binds -- but it is checked and reported.

The de-duplicated topologies are returned as `DemographicTopology` objects
(one canonical event ordering each), ready for the existing Stan pipeline.
"""

from itertools import combinations
import networkx as nx


# ----------------------------------------------------------------------
# State representation
#   nx.DiGraph, edges point ancestor -> descendant (parent -> child).
#   node attrs:  kind  in {'leaf','internal'}
#                label (leaf name, '' for internal)
#                active (bool) -- on the current backward frontier
# ----------------------------------------------------------------------

def _sig(d):
    """Node signature used for isomorphism (kind + leaf label + active)."""
    return f"{d['kind']}|{d.get('label','')}|{int(d['active'])}"


def _wl(g):
    for n, d in g.nodes(data=True):
        d['_sig'] = _sig(d)
    return nx.weisfeiler_lehman_graph_hash(g, node_attr='_sig', iterations=4)


def _iso(g1, g2):
    return nx.is_isomorphic(
        g1, g2, node_match=lambda a, b: _sig(a) == _sig(b))


def _initial(leaves):
    g = nx.DiGraph()
    for x in leaves:
        g.add_node(x, kind='leaf', label=x, active=True)
    return g


def _active(g):
    return [n for n, d in g.nodes(data=True) if d['active']]


def _fresh(g):
    ints = [n for n in g.nodes if isinstance(n, int)]
    return (max(ints) + 1) if ints else 0


def _merge(g, x, y):
    h = g.copy()
    m = _fresh(h)
    h.add_node(m, kind='internal', label='', active=True)
    h.add_edge(m, x); h.add_edge(m, y)
    h.nodes[x]['active'] = False
    h.nodes[y]['active'] = False
    return h


def _admix(g, x):
    h = g.copy()
    p1 = _fresh(h); p2 = p1 + 1
    h.add_node(p1, kind='internal', label='', active=True)
    h.add_node(p2, kind='internal', label='', active=True)
    h.add_edge(p1, x); h.add_edge(p2, x)
    h.nodes[x]['active'] = False
    return h


# ----------------------------------------------------------------------
# Backward generation with iso-class de-duplication of partial states
# ----------------------------------------------------------------------

def _generate(leaves, max_admix):
    start = _initial(leaves)
    # stack of (graph, admix_remaining)
    stack = [(start, max_admix)]
    seen = {}          # (budget, wl) -> [rep graphs]   (avoid re-expanding iso states)
    completed = {}     # wl -> [rep graphs]             (distinct finished topologies)

    def seen_before(store, g, *extra):
        key = (extra, _wl(g))
        bucket = store.setdefault(key, [])
        for r in bucket:
            if _iso(r, g):
                return True
        bucket.append(g)
        return False

    while stack:
        g, budget = stack.pop()
        if seen_before(seen, g, budget):
            continue
        act = _active(g)
        if len(act) == 1:                       # reached the root -> a topology
            seen_before(completed, g)           # de-dup ignoring budget
            continue
        for x, y in combinations(sorted(act, key=str), 2):
            stack.append((_merge(g, x, y), budget))
        if budget > 0:
            for x in act:
                stack.append((_admix(g, x), budget - 1))

    out = []
    for bucket in completed.values():
        out.extend(bucket)
    return out


# ----------------------------------------------------------------------
# Structural classification + pruning
# ----------------------------------------------------------------------

def _roles(g):
    """Return dict node -> role in {'root','leaf','split','admix'} or None."""
    role = {}
    for n, d in g.nodes(data=True):
        ind, outd = g.in_degree(n), g.out_degree(n)
        if d['kind'] == 'leaf':
            role[n] = 'leaf'
        elif ind == 0 and outd == 2:
            role[n] = 'root'
        elif ind == 1 and outd == 2:
            role[n] = 'split'
        elif ind == 2 and outd == 1:
            role[n] = 'admix'
        elif ind == 1 and outd == 1:
            role[n] = 'segment'         # passthrough admixture-source segment
        else:
            role[n] = None
    return role


def _is_trivial_bubble(g):
    """True if some admixture's two sources merge straight back together
    with nothing else attached -> non-identifiable proportion."""
    for a in g.nodes:
        parents = list(g.predecessors(a))
        if len(parents) != 2:
            continue
        p1, p2 = parents
        gp1 = list(g.predecessors(p1))
        gp2 = list(g.predecessors(p2))
        if (g.out_degree(p1) == 1 and g.out_degree(p2) == 1
                and len(gp1) == 1 and gp1 == gp2):
            # p1,p2 are sole-purpose sources sharing one common ancestor m
            return True
    return False


def _n_admix(g):
    return sum(1 for n in g.nodes if g.in_degree(n) == 2)


def _param_count(g, n_bins, model='ibd', nvarying=True):
    """Free parameters of the demographic model for this topology."""
    n_nodes = g.number_of_nodes()
    n_events = sum(1 for n in g.nodes if g.in_degree(n) == 2) \
        + sum(1 for n in g.nodes if g.out_degree(n) == 2 and g.in_degree(n) >= 1)
    # events = #admix + #merge(non-root) ; root merge also an event:
    n_events = (g.number_of_nodes() - len([n for n in g.nodes
                                           if g.nodes[n]['kind'] == 'leaf']))
    n_admix = _n_admix(g)
    ne = n_nodes if nvarying else 1
    times = n_events
    return ne + times + n_admix


def _stat_count(n_leaves, n_bins, model='ibd'):
    if model == 'snp':
        return n_leaves * (n_leaves - 1) // 2
    ibd = n_leaves * (n_leaves + 1) // 2 * n_bins
    if model == 'ibd':
        return ibd
    return ibd + n_leaves * (n_leaves - 1) // 2   # mixed


# ----------------------------------------------------------------------
# Convert a distinct DAG to a DemographicTopology (one canonical ordering)
# ----------------------------------------------------------------------

def _to_demography(g, leaves):
    from demography import DemographicTopology
    dem = DemographicTopology(list(leaves))
    name = {n: n for n in g.nodes if g.nodes[n]['kind'] == 'leaf'}
    counter = [0]

    def newname():
        counter[0] += 1
        return f"I{counter[0]}"

    active = set(n for n, d in g.nodes(data=True) if d['kind'] == 'leaf')
    done = set()

    def emit_admix(child):
        p1, p2 = sorted(g.predecessors(child), key=str)
        for p in (p1, p2):
            name[p] = newname()
        dem.add_admixture_event(name[child], name[p1], name[p2])
        active.discard(child); active.update([p1, p2])
        done.update([p1, p2])

    while len(active) > 1:
        progressed = False
        # prefer a ready MERGE (a node whose both children are active)
        for m in g.nodes:
            if m in done or g.nodes[m]['kind'] == 'leaf':
                continue
            kids = list(g.successors(m))
            if g.out_degree(m) == 2 and all(k in active for k in kids):
                name[m] = newname()
                c1, c2 = sorted(kids, key=lambda k: name[k])
                dem.add_merge_event(name[c1], name[c2], name[m])
                active.difference_update(kids); active.add(m)
                done.add(m)
                progressed = True
                break
        if progressed:
            continue
        # otherwise perform a ready ADMIX (an active node with two parents)
        for x in sorted(active, key=str):
            if g.in_degree(x) == 2 and not all(p in done for p in g.predecessors(x)):
                emit_admix(x)
                progressed = True
                break
        if not progressed:
            raise RuntimeError("could not linearise topology (unexpected structure)")

    dem.finalize_root()
    return dem


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def enumerate_topologies(leaves, max_admix=2, n_bins=10, model='ibd',
                         nvarying=True, drop_bubbles=True, verbose=True):
    """Return a list of (DemographicTopology, info-dict) for every distinct,
    structurally-valid, identifiable-by-counting topology."""
    raw = _generate(list(leaves), max_admix)
    n = len(leaves)
    m = _stat_count(n, n_bins, model)

    kept, by_a = [], {}
    n_bubble = n_overparam = n_badrole = 0
    for g in raw:
        roles = _roles(g)
        if any(r is None for r in roles.values()):
            n_badrole += 1
            continue
        if drop_bubbles and _is_trivial_bubble(g):
            n_bubble += 1
            continue
        p = _param_count(g, n_bins, model, nvarying)
        if p > m:
            n_overparam += 1
            continue
        a = _n_admix(g)
        by_a[a] = by_a.get(a, 0) + 1
        kept.append((g, {'n_admix': a, 'n_params': p, 'n_stats': m}))

    if verbose:
        print(f"leaves={list(leaves)}  max_admix={max_admix}  model={model}")
        print(f"raw finished DAGs (iso-distinct): {len(raw)}")
        print(f"  dropped (bad role / non-binary): {n_badrole}")
        print(f"  dropped (trivial bubble)        : {n_bubble}")
        print(f"  dropped (p>m, p={'?'} m={m})    : {n_overparam}")
        print(f"kept topologies: {len(kept)}")
        for a in sorted(by_a):
            print(f"    {a} admixture event(s): {by_a[a]}")

    return [(_to_demography(g, leaves), info) for g, info in kept]


if __name__ == "__main__":
    import sys
    nleaf = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    names = [chr(ord('A') + i) for i in range(nleaf)]
    res = enumerate_topologies(names, max_admix=kmax, verbose=True)
