# Topology selection

Search over admixture graphs, scored by the per-branch model (`branch_model/`), checked
against a full enumeration on the four-leaf data.

| file | role |
|---|---|
| `graphs.py` | graphs as topologies, canonical keys, the moves, constraints on the space |
| `test_graphs.py` | checks of the space and the moves (`python -B -m unittest test_graphs -v`) |
| `run.py` | `fit` (every graph in the space), `search` (searches on the cached scores), `report` |
| `summarize.py` | `runs/default/report.md` and the figures below |

Figures in `runs/default/`: `landscape.png` (every graph's ELBO relative to the best, per
dataset: the ground truth the searches are measured against), `graphs.png` (the graphs
that matter, drawn with schematic heights, with their ELBO gaps), `search.png` (share of
searches that have found the best graph, by number of fits).

## 1. Graphs and identity

A graph is a set of events, each consuming and creating branches:
`merge (x, y) → p`, `admix x → (y, z)`.  Leaves have no creator; the root has no
consumer.  No heights, no event order.

The key labels every branch by its subnetwork (a leaf by its name, a merge output by its
sorted inputs, an admixture source by `[input]` and a side) and takes the root's label,
minimised over the 2^k side assignments, since an admixture's two sources are
exchangeable (f ↔ 1 − f).  For example, the four-leaf truth is `((([b]1,c),d),([b]2,a))`.

## 2. Moves

Primitives: `detach(b)` (remove the merge consuming b; its partner takes the merge's
place), `attach(b, t)` (new merge of b and t; valid unless t flows into b), `split(x)`
(new admixture on x; one source continues, one dangles), `unsplit`.

| move | built from | reverse |
|---|---|---|
| re-merge | `detach(b)`, `attach(b, t)`: NNI, SPR, re-rooting, moving an admixture source | re-merge b with its old partner |
| add (k < K) | `split(x)`, `attach(z, t)` | remove, dropping z |
| remove | `detach(source)`, `unsplit` | add |
| relocate | remove then add as one move; only the end point is scored | relocate back |

Constraints (checked on every proposed graph, so moves stay reversible): at most K
admixtures; no bubble (the two sources of an admixture merging with each other next);
no ghost (a source joining only above every other lineage); optionally only leaves
admixed.  Always excluded: the root lineage admixed with itself.

Proposal: one move uniformly from all valid moves that change the graph,
`q(G → G′) = m(G, G′) / |M(G)|` with `m` the number of moves landing on G′.
Metropolis–Hastings on cached scores:

    log α = β (ELBO(G′) − ELBO(G)) + log p(G′) − log p(G) + log q(G′ → G) − log q(G → G′)

Prior: uniform over the number of admixtures, then uniform within it.

Tests: tree counts 3 / 15 / 105; in every enumerated space every move is reversible,
every graph satisfies the constraints and is accepted by the per-branch model, and the
moves reach every graph; the key ignores names, listing and source sides.

## 3. Experiment (2026-10-15)

Space: four leaves, at most one admixture, no bubble, no ghost: 153 graphs (15 trees).
Data: T_true subsample 0 of both four-leaf designs at 100 and 1000 cM.  Every graph fitted
with Pathfinder 3 starts × 1 path (153 graphs in 40–50 s with 10 in parallel; no failures).
Searches are replayed on the cached scores, which is exactly what a live search would
fit; "fits" = distinct graphs scored.  Walks: 200 chains × 300 steps from a random tree.
Greedy and stepwise: from each of the 15 trees.

**Best graph.**  At 1000 cM (both designs) and varying 100 cM, the truth is the best graph
by 87–619 nats.  At growth 100 cM the best graph puts b's second source above the (c,d)
merge, not with c; the truth is third, 5.6 nats behind (posterior 0.004).

**Searches** (found = share of runs that scored the best graph; fits = median, 90th percentile):

| search | growth 100 | growth 1000 | varying 100 | varying 1000 |
|---|---|---|---|---|
| walk, β = 1 | 100%, 20 / 36 | 40% | 52% | 42% |
| walk, β = 0.1 | 100% | 37% | 64% | 34% |
| walk + relocate, β = 1 | 100%, 19 / 35 | 100%, 23 / 40 | 100%, 25 / 40 | 100%, 24 / 40 |
| greedy | 100% | 73% | 93% | 67% |
| greedy + relocate | 100% | 100% | 100% | 100% |
| stepwise (best tree, then all) | 100%, 23 | 100%, 21 | 100%, 21 | 100%, 21 |

- **Walks without relocate get trapped** at the best graph with the wrong population
  admixed (a, c, d or (c,d) instead of b), 600–1,800 nats below the truth.  Leaving needs
  a remove (through a tree hundreds of nats worse); tempering to β = 0.1 does not help.
  Relocate removes the barrier: 100% in every dataset, about 25 fits to the best graph.
- **Stepwise** (hill-climb over trees, then over the full space) found the best graph in
  every run with 21–23 fits and stops after ~45–50.
- Here the posterior is on one graph, so any search that scores the best graph covers
  all the mass.  The walks are the route to p(G | x) where graphs are close (growth
  100 cM: 0.99 / 0.006 / 0.004).
