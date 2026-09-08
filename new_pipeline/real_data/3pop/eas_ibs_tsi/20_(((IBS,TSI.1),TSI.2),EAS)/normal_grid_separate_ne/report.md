# `(((IBS,TSI.1),TSI.2),EAS)`

**Normal, separate IBD/SNP Ne, recent grid** | topology 20 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | +258.30 | +- 0.12 (MC) |
| logZ (importance sampling) | +274.92 | |
| ESS of the IS weights | 1.1 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -40.81 | already applied |
| mode kept / MAP start / runtime | 3 / 11 | 25 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 119.4 +- 0.5 | 119.4 |
| 2 | MERGE | TSI.2 + IBS -> n1 | 5.7 +- 0.0 | 125.1 |
| 3 | MERGE | TSI.1 + n1 -> n2 | 214.1 +- 0.3 | 339.3 |
| 4 | MERGE | n2 + EAS -> root | 5.6 +- 0.0 | 344.9 |

## Admixture fraction

**f = 0.750 +- 0.002** (fraction from `TSI.1`; 0.250 from `TSI.2`)

## Recent effective sizes for IBD (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 46,901,978 | 27,632,747 |
| `IBS` | 1,989,707 | 1,162,340 |
| `TSI` | 684,107 | 615,097 |

## Effective sizes for IBD (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 866,803 | 0.01 |
| `IBS` | 523,283 | 0.03 |
| `TSI` | 400,173 | 0.04 |
| `TSI.1` | 49,577 | 0.04 |
| `TSI.2` | 53,816 | 0.02 |
| `n1` | 16,558 | 0.02 |
| `n2` | 97 | 0.02 |
| `root` | 70 | 0.02 |

log-Ne random-walk step scale tau_ibd = 2.521

## Recent effective sizes for SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,315 | 2,294 |
| `IBS` | 312,795 | 313,297 |
| `TSI` | 183,093 | 184,283 |

## Effective sizes for SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,906 | 0.01 |
| `IBS` | 309,650 | 0.01 |
| `TSI` | 195,284 | 0.01 |
| `TSI.1` | 92,375 | 0.02 |
| `TSI.2` | 288,162 | 0.01 |
| `n1` | 275,207 | 0.01 |
| `n2` | 42,350 | 0.00 |
| `root` | 32,585 | 0.00 |

log-Ne random-walk step scale tau_snp = 0.774

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | +557.3 | 222 | 34.15 |
| SNP | +40.1 | 6 | 1.16 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.01 | -0.32 | +0.35 |
| **IBS** | -0.32 | +0.51 | +0.11 |
| **TSI** | +0.35 | +0.11 | -0.77 |

![spectrum](spectrum_fit.png)
