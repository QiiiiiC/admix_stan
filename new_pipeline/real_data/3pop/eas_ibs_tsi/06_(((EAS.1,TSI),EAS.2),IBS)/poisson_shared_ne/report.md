# `(((EAS.1,TSI),EAS.2),IBS)`

**Poisson, shared Ne** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -42821.60 | +- 1.21 (MC) |
| logZ (importance sampling) | -42753.13 | |
| ESS of the IS weights | 3.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 1 | 12 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 276.0 +- 0.6 | 277.0 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.0 +- 0.0 | 278.0 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 279.0 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,228,230 | 0.04 |
| `IBS` | 243,036 | 0.08 |
| `TSI` | 366,372 | 0.17 |
| `EAS.1` | 1,184,872 | 0.04 |
| `EAS.2` | 46 | 0.06 |
| `n1` | 52 | 0.05 |
| `n2` | 57 | 0.03 |
| `root` | 1,081 | 0.02 |

log-Ne random-walk step scale tau = 2.976

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,656.7 | 222 | 3066.07 |
| SNP | -34,833.4 | 6 | 11625.67 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +111.07 | -110.24 | -109.28 |
| **IBS** | -110.24 | +55.13 | +165.94 |
| **TSI** | -109.28 | +165.94 | +52.26 |

![spectrum](spectrum_fit.png)
