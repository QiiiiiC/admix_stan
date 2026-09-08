# `(((EAS,IBS.1),IBS.2),TSI)`

**Poisson, separate IBD/SNP Ne** | topology 10 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3721.15 | +- 0.61 (MC) |
| logZ (importance sampling) | -3681.34 | |
| ESS of the IS weights | 2.2 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 19 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 85.1 +- 3.3 | 85.1 |
| 2 | MERGE | IBS.2 + EAS -> n1 | 48.5 +- 4.1 | 133.6 |
| 3 | MERGE | IBS.1 + n1 -> n2 | 40.1 +- 0.4 | 173.7 |
| 4 | MERGE | n2 + TSI -> root | 1.1 +- 0.0 | 174.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `IBS.1`; 0.000 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,575,371 | 0.06 |
| `IBS` | 849,487 | 0.11 |
| `TSI` | 319,581 | 0.03 |
| `IBS.1` | 42,982 | 0.10 |
| `IBS.2` | 30,389 | 0.05 |
| `n1` | 30,652 | 0.05 |
| `n2` | 97,897 | 0.03 |
| `root` | 85,116 | 0.04 |

log-Ne random-walk step scale tau_ibd = 1.128

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 762 | 0.01 |
| `IBS` | 138,631 | 0.06 |
| `TSI` | 130,096 | 0.07 |
| `IBS.1` | 122,139 | 0.10 |
| `IBS.2` | 8,618 | 0.04 |
| `n1` | 8,000 | 0.04 |
| `n2` | 22,028 | 0.08 |
| `root` | 22,301 | 0.09 |

log-Ne random-walk step scale tau_snp = 0.943

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,554.8 | 222 | 71.62 |
| SNP | +33.8 | 6 | 3.29 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.01 | +0.05 | -0.03 |
| **IBS** | +0.05 | -0.13 | +0.04 |
| **TSI** | -0.03 | +0.04 | +0.01 |

![spectrum](spectrum_fit.png)
