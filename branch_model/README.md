# Per-branch model: the graph as a topology, not an event order

`mixed_model_Nsmooth_poisson` (IBD counts + SNP covariance, smoothed Ne) with the graph
entering as a partial order of events.  Same likelihood and priors as the model in
`methods/stan/`; a different data structure and a different computation.  One topology
is one model: events that are not related are never put in an order.

| file | role |
|---|---|
| `stan/mixed_model_Nsmooth_poisson_branch.stan` | the model |
| `topology.py` | builds the model's data from a `DemographicTopology` |
| `test_branch.py` | equality with the current model (`python -B -m unittest test_branch -v`) |

## 1. Data structure

A graph is a set of **branches** (populations) and **events**.  Each event consumes
branches and creates branches:

- merge: children `x`, `y` → parent `z`
- admixture: child `x` → sources `y` (fraction f) and `z` (fraction 1 − f)

Going back in time, a branch is created by one event (or is a leaf, starting at 0) and
consumed by one event (or is the root).  That is all the model is told:

| per event | meaning |
|---|---|
| `event_admix[e]` | 0 for a merge, else the index of its admixture fraction |
| `event_nodes[e]` | the three branches `x, y, z` above |
| `pred[e]`, `n_pred[e]` | the events that create this event's inputs: it must come after them |
| `open[e]`, `n_open[e]` | the other branches that exist when the event is processed |

| per branch | meaning |
|---|---|
| `start_event[a]` | the event that creates it (0 = leaf) |
| `end_event[a]` | the event that consumes it (0 = root) |
| `ne_parent[a]`, `ne_admix_idx[a]` | where its Ne walk is centred (section 3) |

Events are listed in any order that respects "must come after".  The listing is
bookkeeping only; a different valid listing gives the same answer (tested).

### Example: the four-leaf admixture graph

Branches `a b c d b1 b2 ab cb cbd root` are numbered 1–10.

| e | event | `event_nodes` | `pred` | `open` |
|---|---|---|---|---|
| 1 | b → b1 (f), b2 (1 − f) | 2, 5, 6 | – | a, c, d |
| 2 | a + b1 → ab | 1, 5, 7 | 1 | c, d, b2 |
| 3 | b2 + c → cb | 6, 3, 8 | 1 | d, ab |
| 4 | cb + d → cbd | 8, 4, 9 | 3 | ab |
| 5 | ab + cbd → root | 7, 9, 10 | 2, 4 | – |

`start_event = [0,0,0,0,1,1,2,3,4,5]`, `end_event = [2,1,3,4,2,3,5,4,5,0]`.

Events 2 and {3, 4} are unrelated: neither is in the other's `pred` chain.  The current
model needs three separate models for their three time orders; here it is one.

## 2. Times

`times[e]` is the wait of event e after the latest event it must come after:

    height[e] = max(height[p] for p in pred[e]) + times[e]        (max of nothing = 0)

`times` has the same prior and lower bound as now (`exponential(0.01)`, ≥ 1).  In the
example, heights 20, 60, 100, 200, 400 are waits 20, 40, 80, 100, 200.  Each branch then
has `start = height[start_event]` (0 for a leaf) and `end = height[end_event]`
(`start + T_max` for the root).  Nothing else in the model depends on how unrelated
events compare in time.

For a graph whose events form one chain this is exactly the current prior.  Where events
are parallel it differs: each waits from its own predecessor, not from whichever event
happened to precede it in time.

## 3. Smoothed Ne

The current log-Ne random walk along the flow, unchanged.  The root has
`log Ne = mu_log + sigma_log * Ne_raw`; every other branch is centred on the branch it
flows into (`ne_parent`), and an admixed branch on the f-weighted mean of its two
sources:

    log Ne[a] = log Ne[parent] + tau * sqrt((end[a] − start[a]) / 100) * Ne_raw[a]

The step uses only the branch's own length, so it needs no event order either.

## 4. Likelihood

### 4.1 Each branch on its own

For branch a with start s, end e, size N:

    surv[a]    = exp(−(e − s) / N)                         a pair that enters together leaves uncoalesced
    dens[b, a] = (1 / 50N) ∫_s^e t exp(−(t − s)/N) (exp(−u t/50) − exp(−v t/50)) dt

`dens[b, a]` is the expected number of segments per cM with length in bin b = [u, v) cM
from a pair that enters the branch together.  It has a closed form
`dens[b, a] = g(u) − g(v)`, with `c = u/50 + 1/N` and

    g(u) = [exp(−u s/50) (c s + 1) − surv · exp(−u e/50) (c e + 1)] / (50 N c²)

Every branch starts and ends at an event, so `exp(−u · height/50)` is computed once per
(event, bin edge) and shared.

### 4.2 One sweep over the events

For every pair of branches (x, y) the model keeps an L × L table over leaf pairs:

    T[x, y][i, j] = P(lineage of leaf i is in x, lineage of leaf j is in y, not yet coalesced)

Two lineages can only coalesce in the same branch.  So `T[x, y]` with x ≠ y is a pure
routing probability with no time in it, and time enters only through `surv[x]` when a
branch ends.  This is why no event order is needed.

Start: `T[i, j][i, j] = 1` for all leaves i, j.  Then for each event:

**Merge x, y → z.**  Pairs already together survive their branch; pairs split between x
and y meet now:

    T[z, z] = surv[x] T[x, x] + surv[y] T[y, y] + T[x, y] + T[y, x]
    T[z, w] = T[x, w] + T[y, w],   T[w, z] = T[w, x] + T[w, y]      for every other open branch w

**Admixture x → y (f), z (1 − f).**  With `S = surv[x] T[x, x]`, each lineage chooses a
source independently:

    T[y, y] = f² S,   T[z, z] = (1 − f)² S,   T[y, z] = T[z, y] = f (1 − f) S
    T[y, w] = f T[x, w],   T[z, w] = (1 − f) T[x, w]   (and the mirror images)

After the sweep, `T[a, a]` holds, for every leaf pair, the mass that arrives in branch a
together and uncoalesced.

### 4.3 Expected IBD and the Poisson term

    number[b, (i, j)] = Σ_a dens[b, a] · T[a, a][i, j]

This is one matrix product: (bins × branches) times (branches × leaf pairs).  The
observed count in bin b for pair (i, j) is Poisson with mean
`number × cm × (number of haplotype pairs)`, evaluated in one vectorised call.

### 4.4 SNP term

`route[a, i]` = probability that leaf i's lineage passes through branch a, filled in the
same sweep (merge: add; admixture: multiply by f or 1 − f).  The drift covariance is

    V[i, j] = Σ_{a ≠ root} (end[a] − start[a]) / Ne[a] · route[a, i] route[a, j]

then double-centred and compared with `w_hat` under a normal with sd `w_se`, as now.

## 5. Tests

At identical parameter values the model equals the current one — expected IBD in every
bin and pair, SNP covariance, branch Ne, both log-likelihoods, to ~10 significant
figures — in every time order of: the four-leaf admixture graph (3 orders), the tree
((a,b),(c,d)) with merges at 50/52 both ways, and a two-admixture graph with a loop
(3 orders).  A listing of events that contradicts their times gives the same answer; a
listing that breaks the partial order is rejected.

## 6. Cost and accuracy (measured 2026-10-06)

T_true of both four-leaf designs, 5 subsamples per cell, Pathfinder.

CPU per fit, 3 starts × 1 path × 1000 draws:

| model | 100 cM | 1000 cM |
|---|---:|---:|
| current (`methods/stan`, IBD-fraction integral already dropped) | 2.9 s | 5.4 s |
| per-branch | 1.0 s | 1.9 s |

About 60% of a path is the ELBO estimate Pathfinder makes at every L-BFGS iteration
(25 draws each); the optimisation itself is about a tenth.  The study wrapper
(`methods/fitting.py`) runs 4 paths per start with the same starting values, which
repeats one optimisation four times; relative to that, 1 path per start costs a quarter
and loses a median 0.2–0.6 nats of best-path ELBO (worst 2.2).

Against the current model (3 starts × 4 paths, an earlier version of this model):
best-path ELBO differed by −0.7 to +0.9 nats in the median (−3.5 to +2.9 over 20 pairs),
and posterior means agreed within a few percent for well-identified parameters at
1000 cM, with no systematic shift.
