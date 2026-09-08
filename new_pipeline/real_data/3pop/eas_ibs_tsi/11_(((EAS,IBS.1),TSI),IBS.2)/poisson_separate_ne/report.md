# `(((EAS,IBS.1),TSI),IBS.2)`

**Poisson, separate IBD/SNP Ne** | topology 11 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3687.56 | +- 0.19 (MC) |
| logZ (importance sampling) | -3671.11 | |
| ESS of the IS weights | 7.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 6 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 85.1 +- 0.3 | 85.1 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 47.5 +- 0.2 | 132.6 |
| 3 | MERGE | n1 + TSI -> n2 | 42.2 +- 0.3 | 174.8 |
| 4 | MERGE | IBS.1 + n2 -> root | 1.0 +- 0.0 | 175.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,578,972 | 0.04 |
| `IBS` | 839,484 | 0.05 |
| `TSI` | 320,038 | 0.03 |
| `IBS.1` | 43,082 | 0.05 |
| `IBS.2` | 26,162 | 0.03 |
| `n1` | 32,342 | 0.02 |
| `n2` | 103,027 | 0.01 |
| `root` | 81,976 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.330

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 758 | 0.01 |
| `IBS` | 167,435 | 0.07 |
| `TSI` | 115,478 | 0.18 |
| `IBS.1` | 109,742 | 0.07 |
| `IBS.2` | 8,527 | 0.01 |
| `n1` | 7,904 | 0.01 |
| `n2` | 20,213 | 0.00 |
| `root` | 20,738 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.999

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,539.8 | 222 | 72.17 |
| SNP | +39.6 | 6 | 1.36 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.03 | -0.05 | -0.01 |
| **IBS** | -0.05 | -0.26 | +0.39 |
| **TSI** | -0.01 | +0.39 | -0.35 |

![spectrum](spectrum_fit.png)
