# `(((IBS,TSI.1),EAS),TSI.2)`

**Normal, shared Ne** | topology 19 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -1681.14 | +- 0.36 (MC) |
| logZ (importance sampling) | -1660.33 | |
| ESS of the IS weights | 7.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.2 +- 0.0 | 1.2 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 211.4 +- 0.3 | 212.7 |
| 3 | MERGE | n1 + EAS -> n2 | 190.0 +- 0.3 | 402.7 |
| 4 | MERGE | TSI.1 + n2 -> root | 2,094.6 +- 3.2 | 2,497.2 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 773,110 | 0.02 |
| `IBS` | 382,779 | 0.02 |
| `TSI` | 894,991 | 0.11 |
| `TSI.1` | 1,053,208,089,698 | 0.10 |
| `TSI.2` | 488,008 | 0.10 |
| `n1` | 1,066 | 0.01 |
| `n2` | 4 | 0.05 |
| `root` | 49,275 | 0.00 |

log-Ne random-walk step scale tau = 1.097

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -1,490.8 | 222 | 53.02 |
| SNP | +22.9 | 6 | 6.92 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +1.32 | -1.18 | -1.42 |
| **IBS** | -1.18 | +3.37 | -1.20 |
| **TSI** | -1.42 | -1.20 | +3.85 |

![spectrum](spectrum_fit.png)
