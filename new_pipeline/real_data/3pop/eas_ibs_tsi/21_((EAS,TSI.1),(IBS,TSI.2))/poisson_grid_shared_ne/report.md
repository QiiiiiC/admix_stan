# `((EAS,TSI.1),(IBS,TSI.2))`

**Poisson, shared Ne, recent grid** | topology 21 of 21 | admixed leaf: **TSI** | 4 events, 8 nodes

> Recent Ne grid: 0-5 and 5-10 generations. The first demographic event is explicitly constrained to occur after generation 10 (minimum 11 with the current one-generation gap).

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -882.69 | +- 0.36 (MC) |
| logZ (importance sampling) | -857.68 | |
| ESS of the IS weights | 5.8 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -29.61 | already applied |
| mode kept / MAP start / runtime | 1 / 2 | 31 s |
| mode search | 12/12 MAP starts succeeded | 3 distinct modes evaluated |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | TSI -> TSI.1 + TSI.2 | 87.3 +- 0.8 | 87.3 |
| 2 | MERGE | TSI.1 + EAS -> n1 | 45.0 +- 0.5 | 132.2 |
| 3 | MERGE | TSI.2 + IBS -> n2 | 101.0 +- 1.3 | 233.3 |
| 4 | MERGE | n1 + n2 -> root | 173.0 +- 2.2 | 406.3 |

## Admixture fraction

**f = 0.001 +- 0.000** (fraction from `TSI.1`; 0.999 from `TSI.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Recent effective sizes shared by IBD and SNP (haploid)

| population | 0-5 gen | 5-10 gen |
|---|---:|---:|
| `EAS` | 8,333,666 | 6,143,189 |
| `IBS` | 3,233,284 | 2,191,716 |
| `TSI` | 1,481,533 | 1,053,935 |

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 2,314,538 | 0.04 |
| `IBS` | 295,236 | 0.05 |
| `TSI` | 366,545 | 0.06 |
| `TSI.1` | 50,961 | 0.04 |
| `TSI.2` | 1,302,229 | 0.05 |
| `n1` | 40,085 | 0.04 |
| `n2` | 993 | 0.02 |
| `root` | 108 | 0.01 |

log-Ne random-walk step scale tau = 1.922

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -699.3 | 222 | 39.78 |
| SNP | +29.8 | 6 | 4.61 |

IBD chi2/n uses Poisson counting variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | -0.20 | +0.61 | -0.22 |
| **IBS** | +0.61 | +1.25 | -2.58 |
| **TSI** | -0.22 | -2.58 | +2.84 |

![spectrum](spectrum_fit.png)
