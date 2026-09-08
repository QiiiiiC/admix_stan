# `((EAS,IBS),TSI)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 01 of 21 | tree (no admixture) | 2 events, 5 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -7241.72 | +- 0.46 (MC) |
| logZ (importance sampling) | -7216.18 | |
| ESS of the IS weights | 8.4 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -26.08 | already applied |
| mode kept / MAP start / runtime | 1 / 11 | 17 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | MERGE | EAS + IBS -> n1 | 233.7 +- 0.9 | 233.7 |
| 2 | MERGE | n1 + TSI -> root | 7.1 +- 0.2 | 240.8 |

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 65,419,634 | 29,784,404 |
| `IBS` | 4,562,645 | 3,229,771 |
| `TSI` | 2,343,496 | 1,772,653 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,089,976 | 0.03 |
| `IBS` | 237,082 | 0.03 |
| `TSI` | 314,801 | 0.05 |
| `n1` | 2,109 | 0.04 |
| `root` | 4,603 | 0.02 |

log-Ne random-walk step scale tau_ibd = 2.561

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,100 | 1,207 |
| `IBS` | 181,135 | 184,990 |
| `TSI` | 226,749 | 218,939 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,303 | 0.01 |
| `IBS` | 172,362 | 0.07 |
| `TSI` | 208,077 | 0.06 |
| `n1` | 19,582 | 0.02 |
| `root` | 21,337 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.871

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,033.9 | 222 | 397.80 |
| SNP | +41.4 | 6 | 0.73 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.00 | -0.01 | +0.02 |
| **IBS** | -0.01 | -0.20 | +0.24 |
| **TSI** | +0.02 | +0.24 | -0.26 |

![spectrum](spectrum_fit.png)
