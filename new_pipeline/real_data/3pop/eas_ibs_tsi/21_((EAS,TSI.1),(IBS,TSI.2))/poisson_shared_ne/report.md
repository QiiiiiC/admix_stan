# `((EAS,TSI.1),(IBS,TSI.2))`

**Poisson, shared Ne** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -901.93 | +- 0.33 (MC) |
| logZ (importance sampling) | -884.40 | |
| ESS of the IS weights | 7.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 8 | 19 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 11.0 +- 0.2 | 11.0 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 115.2 +- 0.2 | 126.3 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 104.4 +- 0.9 | 230.7 |
| 4 | MERGE | n1 + n2 -> root | 191.9 +- 0.4 | 422.6 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,585,874 | 0.03 |
| `IBS` | 313,160 | 0.04 |
| `TSI` | 913,367 | 0.02 |
| `TSI.1` | 73,215 | 0.02 |
| `TSI.2` | 399,089 | 0.02 |
| `n1` | 47,008 | 0.03 |
| `n2` | 1,098 | 0.01 |
| `root` | 46 | 0.01 |

log-Ne random-walk step scale tau = 1.376

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -774.6 | 222 | 35.69 |
| SNP | +34.5 | 6 | 3.04 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.21 | +0.51 | -0.10 |
| **IBS** | +0.51 | +1.17 | -2.29 |
| **TSI** | -0.10 | -2.29 | +2.34 |

![spectrum](spectrum_fit.png)
