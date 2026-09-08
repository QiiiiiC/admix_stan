# `(((EAS,TSI.1),TSI.2),IBS)`

**Poisson, shared Ne** | topology 18 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -4967.94 | +- 0.18 (MC) |
| logZ (importance sampling) | -4951.39 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 17 s |
| mode search | 10/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 83.8 +- 0.5 | 83.8 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 207.2 +- 0.5 | 291.0 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 10.2 +- 0.0 | 301.2 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 302.3 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,723,638 | 0.02 |
| `IBS` | 236,990 | 0.04 |
| `TSI` | 447,580 | 0.06 |
| `TSI.1` | 126,639 | 0.04 |
| `TSI.2` | 44 | 0.01 |
| `n1` | 57 | 0.01 |
| `n2` | 1,186 | 0.00 |
| `root` | 823 | 0.00 |

log-Ne random-walk step scale tau = 1.881

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -4,829.4 | 222 | 21318.22 |
| SNP | +35.9 | 6 | 2.58 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.49 | -0.27 | -0.70 |
| **IBS** | -0.27 | -1.78 | +2.46 |
| **TSI** | -0.70 | +2.46 | -0.98 |

![spectrum](spectrum_fit.png)
