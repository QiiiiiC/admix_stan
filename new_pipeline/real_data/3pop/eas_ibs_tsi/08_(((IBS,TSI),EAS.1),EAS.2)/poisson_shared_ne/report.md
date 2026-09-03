# `(((IBS,TSI),EAS.1),EAS.2)`

**Poisson, shared Ne** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6057.35 | +- 0.50 (MC) |
| logZ (importance sampling) | -6023.82 | |
| ESS of the IS weights | 1.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 19 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | IBS + TSI -> n1 | 237.8 +- 0.8 | 238.8 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 8.8 +- 0.2 | 247.6 |
| 4 | MERGE | EAS.1 + n2 -> root | 124.8 +- 3.1 | 372.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,377,912 | 0.03 |
| `IBS` | 308,600 | 0.05 |
| `TSI` | 423,147 | 0.12 |
| `EAS.1` | 898,778 | 0.03 |
| `EAS.2` | 126,214 | 0.03 |
| `n1` | 1,491 | 0.02 |
| `n2` | 718 | 0.03 |
| `root` | 66 | 0.01 |

log-Ne random-walk step scale tau = 0.781

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,799.6 | 222 | 108.50 |
| SNP | +33.0 | 6 | 3.54 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.36 | -0.19 | -0.52 |
| **IBS** | -0.19 | +1.90 | -1.65 |
| **TSI** | -0.52 | -1.65 | +2.53 |

![spectrum](spectrum_fit.png)
