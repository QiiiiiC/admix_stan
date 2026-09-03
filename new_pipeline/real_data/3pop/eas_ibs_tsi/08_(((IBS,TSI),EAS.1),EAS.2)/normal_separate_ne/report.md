# `(((IBS,TSI),EAS.1),EAS.2)`

**Normal, separate IBD/SNP Ne** | topology 08 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +3815.43 | +- 0.37 (MC) |
| logZ (importance sampling) | +3843.90 | |
| ESS of the IS weights | 1.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| seed kept / runtime | 1 | 18 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 127.9 +- 0.7 | 127.9 |
| 2 | MERGE | IBS + TSI -> n1 | 21.3 +- 1.4 | 149.2 |
| 3 | MERGE | EAS.2 + n1 -> n2 | 137.4 +- 6.3 | 286.6 |
| 4 | MERGE | EAS.1 + n2 -> root | 350.2 +- 7.2 | 636.8 |

## Admixture fraction

**f = 0.721 +- 0.001** (fraction from `EAS.1`; 0.279 from `EAS.2`)

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,577,825 | 0.04 |
| `IBS` | 288,341 | 0.05 |
| `TSI` | 421,365 | 0.04 |
| `EAS.1` | 9,838,768 | 0.10 |
| `EAS.2` | 3,414 | 0.03 |
| `n1` | 20,103 | 0.03 |
| `n2` | 3,021 | 0.24 |
| `root` | 18,341 | 0.11 |

log-Ne random-walk step scale tau_ibd = 1.134

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,090 | 0.03 |
| `IBS` | 114,497 | 0.07 |
| `TSI` | 101,743 | 0.06 |
| `EAS.1` | 2,491 | 0.03 |
| `EAS.2` | 14,701 | 0.17 |
| `n1` | 65,698 | 0.07 |
| `n2` | 19,283 | 0.17 |
| `root` | 14,146 | 0.01 |

log-Ne random-walk step scale tau_snp = 0.603

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +3,959.2 | 222 | 3.64 |
| SNP | +34.2 | 6 | 3.15 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.02 | -0.04 | -0.01 |
| **IBS** | -0.04 | -0.04 | +0.12 |
| **TSI** | -0.01 | +0.12 | -0.10 |

![spectrum](spectrum_fit.png)
