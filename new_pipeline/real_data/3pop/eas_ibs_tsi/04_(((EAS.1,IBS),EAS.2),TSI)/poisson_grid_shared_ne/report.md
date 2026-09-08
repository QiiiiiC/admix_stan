# `(((EAS.1,IBS),EAS.2),TSI)`

**Poisson, shared Ne, recent grid** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4813.32 | +- 0.38 (MC) |
| logZ (importance sampling) | -4782.79 | |
| ESS of the IS weights | 1.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 3 / 7 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 292.3 +- 0.5 | 292.3 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 8.4 +- 0.1 | 300.7 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.0 +- 0.0 | 301.7 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 302.7 |

## Admixture fraction

**f = 0.991 +- 0.001** (fraction from `EAS.1`; 0.009 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 20,711,324 | 16,014,461 |
| `IBS` | 6,073,551 | 3,823,891 |
| `TSI` | 2,466,984 | 1,696,697 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,590,444 | 0.02 |
| `IBS` | 225,469 | 0.06 |
| `TSI` | 306,454 | 0.04 |
| `EAS.1` | 52 | 0.01 |
| `EAS.2` | 1,865 | 0.05 |
| `n1` | 1,866 | 0.05 |
| `n2` | 1,697 | 0.02 |
| `root` | 813 | 0.02 |

log-Ne random-walk step scale tau = 2.712

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,615.9 | 222 | 22003.84 |
| SNP | +36.0 | 6 | 2.54 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.40 | -0.21 | -0.58 |
| **IBS** | -0.21 | -0.93 | +1.42 |
| **TSI** | -0.58 | +1.42 | -0.23 |

![spectrum](spectrum_fit.png)
