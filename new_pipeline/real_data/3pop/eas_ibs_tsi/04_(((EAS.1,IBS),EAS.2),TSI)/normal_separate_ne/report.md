# `(((EAS.1,IBS),EAS.2),TSI)`

**Normal, separate IBD/SNP Ne** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +302.61 | +- 0.80 (MC) |
| logZ (importance sampling) | +350.93 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 13 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 121.3 +- 0.5 | 121.3 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.2 +- 0.0 | 122.5 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 23.0 +- 0.6 | 145.5 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 146.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,767,829 | 0.07 |
| `IBS` | 621,909 | 0.12 |
| `TSI` | 330,932 | 0.03 |
| `EAS.1` | 42,561 | 0.02 |
| `EAS.2` | 9,516 | 0.05 |
| `n1` | 9,281 | 0.05 |
| `n2` | 72,697 | 0.01 |
| `root` | 78,905 | 0.01 |

log-Ne random-walk step scale tau_ibd = 0.777

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 670 | 0.03 |
| `IBS` | 126,060 | 0.02 |
| `TSI` | 467,782 | 0.02 |
| `EAS.1` | 168,552 | 0.01 |
| `EAS.2` | 604,100 | 0.01 |
| `n1` | 595,532 | 0.01 |
| `n2` | 357,105 | 0.01 |
| `root` | 544,235 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.444

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +659.0 | 222 | 32.41 |
| SNP | +10.4 | 6 | 11.08 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.37 | +0.76 | -0.03 |
| **IBS** | +0.76 | +0.80 | -2.41 |
| **TSI** | -0.03 | -2.41 | +2.33 |

![spectrum](spectrum_fit.png)
