# `(((EAS.1,IBS),EAS.2),TSI)`

**Normal, shared Ne** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35834.70 | +- 0.05 (MC) |
| logZ (importance sampling) | -35827.84 | |
| ESS of the IS weights | 2.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 3 / 7 | 14 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 189.7 +- 0.2 | 189.7 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.1 +- 0.0 | 190.8 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 2.1 +- 0.0 | 192.9 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 193.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,480,790 | 0.02 |
| `IBS` | 625,772 | 0.03 |
| `TSI` | 328,369 | 0.02 |
| `EAS.1` | 683 | 0.02 |
| `EAS.2` | 95 | 0.02 |
| `n1` | 84 | 0.02 |
| `n2` | 27,472 | 0.01 |
| `root` | 13,782 | 0.01 |

log-Ne random-walk step scale tau = 3.674

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,192.4 | 222 | 49.77 |
| SNP | -34,363.4 | 6 | 11469.01 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +110.39 | -96.79 | -121.48 |
| **IBS** | -96.79 | +48.94 | +145.11 |
| **TSI** | -121.48 | +145.11 | +95.12 |

![spectrum](spectrum_fit.png)
