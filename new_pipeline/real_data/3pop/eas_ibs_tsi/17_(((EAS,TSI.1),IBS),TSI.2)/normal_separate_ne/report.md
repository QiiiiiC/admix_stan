# `(((EAS,TSI.1),IBS),TSI.2)`

**Normal, separate IBD/SNP Ne** | topology 17 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -493.14 | +- 0.92 (MC) |
| logZ (importance sampling) | -442.66 | |
| ESS of the IS weights | 4.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 1 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 136.5 +- 0.2 | 136.5 |
| 2 | MERGE | TSI.2 + EAS -> n1 | 1.9 +- 0.0 | 138.5 |
| 3 | MERGE | n1 + IBS -> n2 | 9.7 +- 0.1 | 148.1 |
| 4 | MERGE | TSI.1 + n2 -> root | 1.0 +- 0.0 | 149.1 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `TSI.1`; 0.000 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,589,077 | 0.02 |
| `IBS` | 218,323 | 0.02 |
| `TSI` | 362,388 | 0.06 |
| `TSI.1` | 53,237 | 0.02 |
| `TSI.2` | 12,814 | 0.03 |
| `n1` | 12,814 | 0.03 |
| `n2` | 75,308 | 0.01 |
| `root` | 70,102 | 0.01 |

log-Ne random-walk step scale tau_ibd = 0.540

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 761 | 0.03 |
| `IBS` | 644,373,097,910 | 0.02 |
| `TSI` | 86,148,227,457 | 0.01 |
| `TSI.1` | 3,633,159,370 | 0.01 |
| `TSI.2` | 81,794,006 | 0.01 |
| `n1` | 81,794,706 | 0.01 |
| `n2` | 527,730,242 | 0.01 |
| `root` | 1,518,088,799 | 0.00 |

log-Ne random-walk step scale tau_snp = 1.363

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +201.0 | 222 | 36.43 |
| SNP | -13.8 | 6 | 19.14 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.71 | +0.76 | +0.64 |
| **IBS** | +0.76 | +2.94 | -4.70 |
| **TSI** | +0.64 | -4.70 | +3.20 |

![spectrum](spectrum_fit.png)
