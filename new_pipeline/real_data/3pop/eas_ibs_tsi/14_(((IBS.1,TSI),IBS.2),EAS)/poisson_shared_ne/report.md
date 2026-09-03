# `(((IBS.1,TSI),IBS.2),EAS)`

**Poisson, shared Ne** | topology 14 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6212.03 | +- 0.70 (MC) |
| logZ (importance sampling) | -6169.54 | |
| ESS of the IS weights | 4.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 242.0 +- 0.9 | 242.0 |
| 2 | MERGE | IBS.2 + TSI -> n1 | 1.0 +- 0.0 | 243.0 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 1.0 +- 0.0 | 244.0 |
| 4 | MERGE | n2 + EAS -> root | 119.6 +- 4.3 | 363.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 887,499 | 0.02 |
| `IBS` | 326,920 | 0.06 |
| `TSI` | 428,786 | 0.11 |
| `IBS.1` | 742 | 0.04 |
| `IBS.2` | 853 | 0.04 |
| `n1` | 846 | 0.04 |
| `n2` | 671 | 0.04 |
| `root` | 123 | 0.24 |

log-Ne random-walk step scale tau = 0.614

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,795.1 | 222 | 113.40 |
| SNP | +11.3 | 6 | 10.78 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.80 | +0.08 | -1.67 |
| **IBS** | +0.08 | -4.32 | +4.46 |
| **TSI** | -1.67 | +4.46 | -1.01 |

![spectrum](spectrum_fit.png)
