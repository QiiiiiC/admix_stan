# `(((EAS,IBS.1),IBS.2),TSI)`

**Poisson, shared Ne, recent grid** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4823.27 | +- 0.30 (MC) |
| logZ (importance sampling) | -4804.57 | |
| ESS of the IS weights | 7.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| seed kept / runtime | 7 | 20 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 291.3 +- 0.6 | 291.3 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 1.0 +- 0.0 | 292.3 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 9.6 +- 0.0 | 302.0 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 303.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 22,753,039 | 13,534,224 |
| `IBS` | 3,432,773 | 2,136,226 |
| `TSI` | 2,910,006 | 1,850,831 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,586,049 | 0.01 |
| `IBS` | 225,437 | 0.03 |
| `TSI` | 307,319 | 0.02 |
| `IBS.1` | 16,847 | 0.01 |
| `IBS.2` | 54 | 0.01 |
| `n1` | 54 | 0.01 |
| `n2` | 1,494 | 0.01 |
| `root` | 798 | 0.01 |

log-Ne random-walk step scale tau = 2.586

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,616.1 | 222 | 22158.13 |
| SNP | +31.8 | 6 | 3.94 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.31 | -0.18 | -0.43 |
| **IBS** | -0.18 | -1.05 | +1.48 |
| **TSI** | -0.43 | +1.48 | -0.57 |

![spectrum](spectrum_fit.png)
