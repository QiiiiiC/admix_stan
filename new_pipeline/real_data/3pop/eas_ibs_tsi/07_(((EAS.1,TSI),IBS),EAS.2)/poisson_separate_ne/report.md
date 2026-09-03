# `(((EAS.1,TSI),IBS),EAS.2)`

**Poisson, separate IBD/SNP Ne** | topology 07 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5900.41 | +- 0.72 (MC) |
| logZ (importance sampling) | -5851.08 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 158.1 +- 3.0 | 159.1 |
| 3 | MERGE | n1 + IBS -> n2 | 1.0 +- 0.0 | 160.1 |
| 4 | MERGE | EAS.1 + n2 -> root | 160.4 +- 3.5 | 320.5 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 901,072 | 0.03 |
| `IBS` | 314,398 | 0.04 |
| `TSI` | 426,484 | 0.05 |
| `EAS.1` | 898,143 | 0.03 |
| `EAS.2` | 15,887 | 0.11 |
| `n1` | 16,786 | 0.12 |
| `n2` | 14,320 | 0.10 |
| `root` | 606 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.291

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 833,955 | 0.72 |
| `IBS` | 116,281 | 0.07 |
| `TSI` | 204,078 | 0.05 |
| `EAS.1` | 847,240 | 0.72 |
| `EAS.2` | 927 | 0.02 |
| `n1` | 1,083 | 0.04 |
| `n2` | 892 | 0.03 |
| `root` | 10,434 | 0.19 |

log-Ne random-walk step scale tau_snp = 1.340

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,682.7 | 222 | 46.45 |
| SNP | +23.7 | 6 | 6.63 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.02 | -0.07 | +0.12 |
| **IBS** | -0.07 | -0.37 | +0.54 |
| **TSI** | +0.12 | +0.54 | -0.74 |

![spectrum](spectrum_fit.png)
