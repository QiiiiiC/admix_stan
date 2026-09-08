# `(((EAS.1,IBS),TSI),EAS.2)`

**Normal, shared Ne, recent grid** | topology 05 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3994.52 | +- 0.08 (MC) |
| logZ (importance sampling) | +4003.67 | |
| ESS of the IS weights | 7.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 2 / 7 | 39 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 130.0 +- 0.4 | 130.0 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 4.0 +- 0.0 | 134.0 |
| 3 | MERGE | n1 + TSI -> n2 | 19.9 +- 0.3 | 154.0 |
| 4 | MERGE | EAS.1 + n2 -> root | 2,311.6 +- 6.3 | 2,465.6 |

## Admixture fraction

**f = 0.992 +- 0.000** (fraction from `EAS.1`; 0.008 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 6,665,628 | 5,451,659 |
| `IBS` | 2,993,299 | 2,036,788 |
| `TSI` | 681,492 | 577,692 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,373,803 | 0.04 |
| `IBS` | 493,438 | 0.03 |
| `TSI` | 401,933 | 0.02 |
| `EAS.1` | 41,481 | 0.01 |
| `EAS.2` | 8,922 | 0.01 |
| `n1` | 9,245 | 0.01 |
| `n2` | 18,228 | 0.01 |
| `root` | 63,782 | 0.00 |

log-Ne random-walk step scale tau = 1.835

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +4,156.5 | 222 | 1.38 |
| SNP | +37.8 | 6 | 1.94 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | +0.99 | -1.06 |
| **IBS** | +0.99 | -2.02 | +0.14 |
| **TSI** | -1.06 | +0.14 | +1.90 |

![spectrum](spectrum_fit.png)
