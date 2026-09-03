# `(((IBS,TSI),EAS.1),EAS.2)`

**Poisson, separate IBD/SNP Ne** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -5860.61 | +- 0.43 (MC) |
| logZ (importance sampling) | -5831.98 | |
| ESS of the IS weights | 1.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 7 | 17 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | IBS + TSI -> n1 | 160.0 +- 1.4 | 161.0 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 26.5 +- 0.2 | 187.5 |
| 4 | MERGE | EAS.1 + n2 -> root | 142.2 +- 0.6 | 329.7 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,112,097 | 0.04 |
| `IBS` | 314,973 | 0.05 |
| `TSI` | 423,941 | 0.02 |
| `EAS.1` | 899,677 | 0.04 |
| `EAS.2` | 98,871 | 0.06 |
| `n1` | 13,960 | 0.02 |
| `n2` | 13,633 | 0.02 |
| `root` | 418 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.509

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 298,288 | 0.09 |
| `IBS` | 137,286 | 0.17 |
| `TSI` | 98,369 | 0.14 |
| `EAS.1` | 253,864 | 0.09 |
| `EAS.2` | 4,601 | 0.05 |
| `n1` | 2,339 | 0.01 |
| `n2` | 848 | 0.01 |
| `root` | 9,255 | 0.00 |

log-Ne random-walk step scale tau_snp = 1.289

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -5,683.4 | 222 | 46.46 |
| SNP | +38.4 | 6 | 1.74 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.00 | -0.16 | +0.16 |
| **IBS** | -0.16 | +0.09 | +0.24 |
| **TSI** | +0.16 | +0.24 | -0.52 |

![spectrum](spectrum_fit.png)
