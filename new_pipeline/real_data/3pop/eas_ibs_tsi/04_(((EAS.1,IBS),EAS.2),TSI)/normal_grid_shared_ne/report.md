# `(((EAS.1,IBS),EAS.2),TSI)`

**Normal, shared Ne, recent grid** | topology 04 of 21 | admixed leaf: **EAS** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -35359.76 | +- 0.23 (MC) |
| logZ (importance sampling) | -35339.83 | |
| ESS of the IS weights | 1.6 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 7 | 40 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | EAS -> EAS.1 + EAS.2 | 244.3 +- 0.3 | 244.3 |
| 2 | MERGE | EAS.2 + IBS -> n1 | 1.0 +- 0.0 | 245.3 |
| 3 | MERGE | EAS.1 + n1 -> n2 | 1.0 +- 0.0 | 246.3 |
| 4 | MERGE | n2 + TSI -> root | 1.0 +- 0.0 | 247.3 |

## Admixture fraction

**f = 0.999 +- 0.000** (fraction from `EAS.1`; 0.001 from `EAS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 2,515,169 | 905,590 |
| `IBS` | 44,091,746 | 14,986,999 |
| `TSI` | 4,445,337 | 3,031,023 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 3,842,059 | 0.05 |
| `IBS` | 283,810 | 0.05 |
| `TSI` | 365,020 | 0.04 |
| `EAS.1` | 20 | 0.02 |
| `EAS.2` | 77 | 0.03 |
| `n1` | 75 | 0.03 |
| `n2` | 5,235 | 0.03 |
| `root` | 780 | 0.02 |

log-Ne random-walk step scale tau = 4.207

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -27,910.4 | 222 | 290.84 |
| SNP | -7,031.1 | 6 | 2358.25 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +49.72 | -42.46 | -55.86 |
| **IBS** | -42.46 | +17.20 | +68.23 |
| **TSI** | -55.86 | +68.23 | +42.32 |

![spectrum](spectrum_fit.png)
