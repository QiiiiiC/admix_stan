# Relaxed graph: the graph as the only random object

A prototype for topology search without fitting every graph separately.  The
parameters (event times, admixture fraction, one haploid Ne) are **unknown
constants**, estimated by maximum likelihood; the **graph** is the random
object.  The full mathematics is in [MATHS.md](MATHS.md).

- **Skeleton**: [admixture, merge, ..., merge, root] on L leaves.  The last slot is a
  compulsory root that moves every deme into one, so every pair coalesces.  For 3
  leaves: 108 slot choices ("vertices") = 45 distinct graphs.
- **Exact EM** (`stan/discrete.stan`): maximise log Σ_v p(v) exp ℓ(θ, v) over θ; the graph
  posterior at θ̂ is the softmax over vertices.  Costs one likelihood per vertex.
- **Profile** (`discrete.stan` with one vertex): each graph with its own θ̂.
- **Relaxed EM** (`stan/relaxed_estep.stan`, `stan/relaxed_mstep.stan`): Gumbel(0,1) logits per
  slot, soft slot matrices M_k = Σ_c softmax(y_k/λ)_c M_k^(c); Monte Carlo EM with NUTS on
  the logits (θ as data) and L-BFGS on θ.  As λ → 0 it tends to the exact EM.  Costs one
  likelihood per evaluation whatever the number of graphs.

All models share `stan/graph.stanfunctions`; at every binary vertex its likelihood equals
`methods/stan/mixed_Nfixed_poisson.stan` on the same graph (tested).

## Data

Borrowed: the tree sequences of `new_pipeline/Nfixed/growth_asymmetric_loop`
(`constant` scenario: 3 leaves, haploid Ne 15,000 on every branch, 50 × 50 cM blocks,
15 diploids per leaf).  IBD and SNP sums are re-extracted in 22 bins from 1 cM.  The
generating graph has a loop on b (10–25 gen) and an admixture of b (30 gen) — two
admixtures — so it is **outside** the one-admixture skeleton; the closest graph is
"admix b; a+b1; b2+c; root".

## Run

```sh
python -B run.py extract                       # block caches from the borrowed trees
python -B run.py discrete                      # exact EM, every genome length in config.json
python -B run.py profile
python -B run.py relaxed --lambdas 1 0.5       # relaxed Monte Carlo EM per temperature
python -B run.py summarize                     # runs/default/summary/{report.md, *.png}
python -B -m unittest test_relaxed -v
```
