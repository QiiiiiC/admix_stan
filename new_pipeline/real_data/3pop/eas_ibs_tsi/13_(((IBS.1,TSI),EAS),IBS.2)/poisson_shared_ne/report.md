# `(((IBS.1,TSI),EAS),IBS.2)`

**Poisson, shared Ne** | topology 13 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5593.17 | +- 0.30 (MC) |
| logZ (importance sampling) | -5573.76 | |
| ESS of the IS weights | 6.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 8 | 19 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 58.3 +- 0.7 | 58.3 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 183.2 +- 0.6 | 241.5 |
| 3 | MERGE | n1 + EAS -> n2 | 157.0 +- 0.5 | 398.5 |
| 4 | MERGE | IBS.1 + n2 -> root | 218.6 +- 0.6 | 617.1 |

## Admixture fraction

**f = 0.003 +- 0.000** (fraction from `IBS.1`; 0.997 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 896,562 | 0.02 |
| `IBS` | 1,054,762 | 0.06 |
| `TSI` | 406,862 | 0.02 |
| `IBS.1` | 4,793 | 0.01 |
| `IBS.2` | 111,001 | 0.04 |
| `n1` | 870 | 0.02 |
| `n2` | 8 | 0.02 |
| `root` | 1,727 | 0.01 |

log-Ne random-walk step scale tau = 1.788

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,475.0 | 222 | 112.25 |
| SNP | +35.7 | 6 | 2.63 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.01 | +0.23 | -0.25 |
| **IBS** | +0.23 | -0.10 | -0.35 |
| **TSI** | -0.25 | -0.35 | +0.80 |

![spectrum](spectrum_fit.png)
