# `(((EAS.1,IBS),EAS.2),TSI)`

**Poisson, separate IBD/SNP Ne, recent grid** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -6985.24 | +- 1.70 (MC) |
| logZ (importance sampling) | -6870.22 | |
| ESS of the IS weights | 3.5 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| seed kept / runtime | 1 | 22 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 11.0 +- 0.0 | 11.0 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 186.4 +- 4.7 | 197.4 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 28.5 +- 0.9 | 226.0 |
| 4 | MERGE | n2 + TSI -> root | 4.9 +- 0.4 | 230.9 |

## Admixture fraction

**f = 1.000 +- 0.000** (fraction from `EAS.1`; 0.000 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 40,839,240 | 24,913,031 |
| `IBS` | 2,651,971 | 1,728,471 |
| `TSI` | 1,755,154 | 1,239,122 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 3,368,412 | 0.08 |
| `IBS` | 418,457 | 0.10 |
| `TSI` | 318,444 | 0.05 |
| `EAS.1` | 1,050,790 | 0.06 |
| `EAS.2` | 71 | 0.17 |
| `n1` | 1,162 | 0.15 |
| `n2` | 2,857 | 0.10 |
| `root` | 6,477 | 0.16 |

log-Ne random-walk step scale tau_ibd = 2.037

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 1,856 | 1,704 |
| `IBS` | 335,416 | 307,560 |
| `TSI` | 312,039 | 296,904 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,354 | 0.03 |
| `IBS` | 243,111 | 0.04 |
| `TSI` | 251,208 | 0.05 |
| `EAS.1` | 1,235 | 0.03 |
| `EAS.2` | 19,219 | 0.03 |
| `n1` | 28,851 | 0.03 |
| `n2` | 32,573 | 0.03 |
| `root` | 37,035 | 0.02 |

log-Ne random-walk step scale tau_snp = 0.475

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -6,673.7 | 222 | 258.54 |
| SNP | +21.7 | 6 | 7.32 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.12 | +0.53 | -0.30 |
| **IBS** | +0.53 | -1.09 | +0.09 |
| **TSI** | -0.30 | +0.09 | +0.48 |

![spectrum](spectrum_fit.png)
