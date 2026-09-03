# `(((IBS,TSI.1),TSI.2),EAS)`

**Poisson, shared Ne** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6113.64 | +- 0.30 (MC) |
| logZ (importance sampling) | -6094.20 | |
| ESS of the IS weights | 11.3 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 7 | 16 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 248.6 +- 1.1 | 249.6 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 1.0 +- 0.0 | 250.6 |
| 4 | MERGE | n2 + EAS -> root | 95.7 +- 2.0 | 346.4 |

## Admixture fraction

**f = 0.000 +- 0.000** (fraction from `TSI.1`; 1.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 892,433 | 0.03 |
| `IBS` | 312,065 | 0.03 |
| `TSI` | 464,615 | 0.04 |
| `TSI.1` | 1,151 | 0.02 |
| `TSI.2` | 424,704 | 0.04 |
| `n1` | 757 | 0.02 |
| `n2` | 534 | 0.02 |
| `root` | 234 | 0.01 |

log-Ne random-walk step scale tau = 0.682

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,805.4 | 222 | 140.43 |
| SNP | +35.9 | 6 | 2.57 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.13 | +0.29 | -0.03 |
| **IBS** | +0.29 | +1.38 | -2.08 |
| **TSI** | -0.03 | -2.08 | +2.02 |

![spectrum](spectrum_fit.png)
