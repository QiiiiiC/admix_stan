# `(((EAS,IBS.1),IBS.2),TSI)`

**Normal, separate IBD/SNP Ne** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +550.75 | +- 0.11 (MC) |
| logZ (importance sampling) | +562.08 | |
| ESS of the IS weights | 3.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 20 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 86.5 +- 0.5 | 86.5 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 43.1 +- 0.3 | 129.6 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 15.2 +- 0.1 | 144.8 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 145.9 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `IBS.1`; 0.001 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,678,449 | 0.03 |
| `IBS` | 829,065 | 0.05 |
| `TSI` | 328,962 | 0.02 |
| `IBS.1` | 39,373 | 0.02 |
| `IBS.2` | 39,361 | 0.02 |
| `n1` | 23,729 | 0.02 |
| `n2` | 103,219 | 0.01 |
| `root` | 80,536 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.420

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 726 | 0.01 |
| `IBS` | 166,903 | 0.06 |
| `TSI` | 123,532 | 0.10 |
| `IBS.1` | 69,523 | 0.05 |
| `IBS.2` | 3,619 | 0.02 |
| `n1` | 8,105 | 0.00 |
| `n2` | 16,335 | 0.00 |
| `root` | 17,005 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.983

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +700.1 | 222 | 32.07 |
| SNP | +41.3 | 6 | 0.76 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.05 | +0.16 | -0.06 |
| **IBS** | +0.16 | -0.13 | -0.20 |
| **TSI** | -0.06 | -0.20 | +0.30 |

![spectrum](spectrum_fit.png)
