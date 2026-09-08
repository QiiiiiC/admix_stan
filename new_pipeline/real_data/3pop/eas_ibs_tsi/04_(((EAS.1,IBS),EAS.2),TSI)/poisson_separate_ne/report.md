# `(((EAS.1,IBS),EAS.2),TSI)`

**Poisson, separate IBD/SNP Ne** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -3733.26 | +- 0.47 (MC) |
| logZ (importance sampling) | -3699.29 | |
| ESS of the IS weights | 1.0 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.78 | already applied |
| mode kept / MAP start / runtime | 1 / 5 | 18 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 126.9 +- 0.6 | 126.9 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 2.3 +- 0.1 | 129.2 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 43.5 +- 0.5 | 172.7 |
| 4 | MERGE | n2 + TSI -> root | 1.1 +- 0.0 | 173.8 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,648,576 | 0.05 |
| `IBS` | 608,214 | 0.09 |
| `TSI` | 319,437 | 0.04 |
| `EAS.1` | 39,888 | 0.02 |
| `EAS.2` | 10,323 | 0.04 |
| `n1` | 10,255 | 0.04 |
| `n2` | 130,218 | 0.01 |
| `root` | 87,822 | 0.02 |

log-Ne random-walk step scale tau_ibd = 1.661

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 773 | 0.01 |
| `IBS` | 219,410 | 0.04 |
| `TSI` | 120,783 | 0.10 |
| `EAS.1` | 2,846 | 0.01 |
| `EAS.2` | 50,866 | 0.02 |
| `n1` | 52,945 | 0.02 |
| `n2` | 17,481 | 0.01 |
| `root` | 15,538 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.931

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -3,567.1 | 222 | 70.61 |
| SNP | +35.8 | 6 | 2.61 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +0.04 | -0.02 | -0.05 |
| **IBS** | -0.02 | -0.28 | +0.35 |
| **TSI** | -0.05 | +0.35 | -0.24 |

![spectrum](spectrum_fit.png)
