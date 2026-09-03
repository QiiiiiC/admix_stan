# `(((IBS,TSI),EAS.1),EAS.2)`

**Normal, shared Ne** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3869.03 | +- 0.20 (MC) |
| logZ (importance sampling) | +3882.00 | |
| ESS of the IS weights | 3.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 126.4 +- 0.5 | 126.4 |
| 2 | MERGE | IBS + TSI -> n1 | 24.4 +- 0.7 | 150.8 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 154.0 +- 0.8 | 304.8 |
| 4 | MERGE | EAS.1 + n2 -> root | 298.4 +- 0.7 | 603.3 |

## Admixture fraction

**f = 0.785 +- 0.001** (fraction from `EAS.1`; 0.215 from `EAS.2`)

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,591,631 | 0.03 |
| `IBS` | 288,757 | 0.03 |
| `TSI` | 420,007 | 0.04 |
| `EAS.1` | 43,638 | 0.02 |
| `EAS.2` | 6,172 | 0.01 |
| `n1` | 19,280 | 0.01 |
| `n2` | 1,112 | 0.01 |
| `root` | 12,053 | 0.01 |

log-Ne random-walk step scale tau = 1.250

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,965.1 | 222 | 3.64 |
| SNP | +29.4 | 6 | 4.75 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.37 | +0.51 | +0.23 |
| **IBS** | +0.51 | +1.89 | -3.07 |
| **TSI** | +0.23 | -3.07 | +2.45 |

![spectrum](spectrum_fit.png)
