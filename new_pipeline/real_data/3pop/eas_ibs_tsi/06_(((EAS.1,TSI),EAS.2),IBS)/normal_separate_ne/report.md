# `(((EAS.1,TSI),EAS.2),IBS)`

**Normal, separate IBD/SNP Ne** | topology 06 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -27.33 | +- 0.23 (MC) |
| logZ (importance sampling) | -9.29 | |
| ESS of the IS weights | 10.9 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 130.6 +- 0.3 | 130.6 |
| 2 | MERGE | EAS.2 + TSI -> n1 | 14.4 +- 0.0 | 145.0 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.0 +- 0.0 | 146.0 |
| 4 | MERGE | n2 + IBS -> root | 1.0 +- 0.0 | 147.1 |

## Admixture fraction

**f = 0.003 +- 0.000** (fraction from `EAS.1`; 0.997 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,676,271 | 0.04 |
| `IBS` | 219,548 | 0.03 |
| `TSI` | 331,071 | 0.02 |
| `EAS.1` | 531,986 | 0.02 |
| `EAS.2` | 22,212 | 0.02 |
| `n1` | 183,350 | 0.01 |
| `n2` | 138,832 | 0.01 |
| `root` | 76,099 | 0.01 |

log-Ne random-walk step scale tau_ibd = 1.786

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,064 | 0.02 |
| `IBS` | 119,914 | 0.08 |
| `TSI` | 208,978 | 0.05 |
| `EAS.1` | 10,790 | 0.01 |
| `EAS.2` | 122 | 0.01 |
| `n1` | 3,031 | 0.01 |
| `n2` | 3,412 | 0.01 |
| `root` | 3,870 | 0.01 |

log-Ne random-walk step scale tau_snp = 1.530

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +182.2 | 222 | 36.55 |
| SNP | +34.4 | 6 | 3.06 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.04 | +0.67 | -0.60 |
| **IBS** | +0.67 | -1.04 | -0.24 |
| **TSI** | -0.60 | -0.24 | +1.37 |

![spectrum](spectrum_fit.png)
