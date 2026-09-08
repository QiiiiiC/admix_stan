# `((EAS,TSI.1),(IBS,TSI.2))`

**Poisson, separate IBD/SNP Ne** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -829.99 | +- 0.29 (MC) |
| logZ (importance sampling) | -810.24 | |
| ESS of the IS weights | 6.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 2 / 3 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 49.9 +- 1.8 | 49.9 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 77.1 +- 1.9 | 126.9 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 37.0 +- 0.8 | 163.9 |
| 4 | MERGE | n1 + n2 -> root | 111.7 +- 3.3 | 275.6 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,579,739 | 0.04 |
| `IBS` | 315,707 | 0.02 |
| `TSI` | 411,585 | 0.06 |
| `TSI.1` | 37,406 | 0.15 |
| `TSI.2` | 456,463 | 0.05 |
| `n1` | 45,851 | 0.02 |
| `n2` | 11,574 | 0.03 |
| `root` | 16,598 | 0.09 |

log-Ne random-walk step scale tau_ibd = 1.296

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,023 | 0.02 |
| `IBS` | 125,266 | 0.05 |
| `TSI` | 127,146 | 0.09 |
| `TSI.1` | 2,320 | 0.17 |
| `TSI.2` | 103,861 | 0.04 |
| `n1` | 2,775 | 0.04 |
| `n2` | 43,058 | 0.20 |
| `root` | 11,237 | 0.35 |

log-Ne random-walk step scale tau_snp = 0.827

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -688.4 | 222 | 3.99 |
| SNP | +39.4 | 6 | 1.40 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.01 | +0.01 | -0.02 |
| **IBS** | +0.01 | -0.15 | +0.14 |
| **TSI** | -0.02 | +0.14 | -0.09 |

![spectrum](spectrum_fit.png)
