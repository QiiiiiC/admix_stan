# `((EAS.1,IBS),(EAS.2,TSI))`

**Poisson, shared Ne** | topology 09 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5049.92 | +- 0.33 (MC) |
| logZ (importance sampling) | -5030.30 | |
| ESS of the IS weights | 3.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 2 / 0 | 17 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 291.2 +- 0.4 | 291.2 |
| 2 | MERGE | EAS.1 + IBS -> n1 | 8.5 +- 0.1 | 299.7 |
| 3 | MERGE | EAS.2 + TSI -> n2 | 1.3 +- 0.0 | 301.0 |
| 4 | MERGE | n1 + n2 -> root | 1.0 +- 0.0 | 302.0 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `EAS.1`; 0.999 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,724,852 | 0.03 |
| `IBS` | 237,495 | 0.02 |
| `TSI` | 323,216 | 0.08 |
| `EAS.1` | 1,809 | 0.03 |
| `EAS.2` | 54 | 0.01 |
| `n1` | 2,268 | 0.00 |
| `n2` | 2,002 | 0.00 |
| `root` | 821 | 0.00 |

log-Ne random-walk step scale tau = 2.503

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,893.8 | 222 | 20739.78 |
| SNP | +32.5 | 6 | 3.72 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.41 | +0.55 | -1.36 |
| **IBS** | +0.55 | -2.82 | +1.90 |
| **TSI** | -1.36 | +1.90 | +0.81 |

![spectrum](spectrum_fit.png)
