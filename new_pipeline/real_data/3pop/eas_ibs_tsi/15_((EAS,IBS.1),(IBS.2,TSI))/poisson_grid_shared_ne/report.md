# `((EAS,IBS.1),(IBS.2,TSI))`

**Poisson, shared Ne, recent grid** | topology 15 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -737.49 | +- 0.71 (MC) |
| logZ (importance sampling) | -697.15 | |
| ESS of the IS weights | 6.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 3 | 27 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 52.0 +- 0.3 | 52.0 |
| 2 | MERGE | IBS.1 + EAS -> n1 | 81.4 +- 1.1 | 133.4 |
| 3 | MERGE | IBS.2 + TSI -> n2 | 109.9 +- 1.8 | 243.3 |
| 4 | MERGE | n1 + n2 -> root | 133.5 +- 2.4 | 376.8 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `IBS.1`; 1.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 13,270,246 | 12,351,729 |
| `IBS` | 2,941,843 | 2,346,422 |
| `TSI` | 9,768,153 | 3,915,794 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,199,133 | 0.03 |
| `IBS` | 1,148,747 | 0.12 |
| `TSI` | 388,581 | 0.09 |
| `IBS.1` | 38,993 | 0.05 |
| `IBS.2` | 129,374 | 0.07 |
| `n1` | 39,453 | 0.06 |
| `n2` | 766 | 0.02 |
| `root` | 379 | 0.00 |

log-Ne random-walk step scale tau = 3.648

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -513.3 | 222 | 73.89 |
| SNP | +36.9 | 6 | 2.24 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.09 | +0.41 | -0.60 |
| **IBS** | +0.41 | -0.08 | -0.76 |
| **TSI** | -0.60 | -0.76 | +1.85 |

![spectrum](spectrum_fit.png)
