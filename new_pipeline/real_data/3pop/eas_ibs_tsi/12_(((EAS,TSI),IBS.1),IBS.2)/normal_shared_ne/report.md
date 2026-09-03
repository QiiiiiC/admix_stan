# `(((EAS,TSI),IBS.1),IBS.2)`

**Normal, shared Ne** | topology 12 of 21 | admixed leaf: **IBS** | 4 events, 8 nodes

> Graph where the two NON-admixed leaves merge first. Excluded by the 'first merge must involve an admixture branch' rule; included here because it is the standard local-clade + deep-source scenario.

## Model score

| quantity | value | |
|---|---|---|
| ELBO | -47491.28 | +- 0.92 (MC) |
| logZ (importance sampling) | -47428.56 | |
| ESS of the IS weights | 1.7 / 4000 | LOW -- treat logZ with suspicion |
| Stan dropped-constant correction | -24.10 | already applied |
| seed kept / runtime | 13 | 6 s |

## Events (in temporal order, most recent first)

| # | type | detail | time (gen) | cumulative |
|---|---|---|---|---|
| 1 | ADMIXTURE | IBS -> IBS.1 + IBS.2 | 1.0 +- 0.0 | 1.0 |
| 2 | MERGE | EAS + TSI -> n1 | 1.0 +- 0.0 | 2.0 |
| 3 | MERGE | IBS.2 + n1 -> n2 | 1.0 +- 0.0 | 3.0 |
| 4 | MERGE | IBS.1 + n2 -> root | 181.9 +- 8.1 | 184.9 |

## Admixture fraction

**f = 0.998 +- 0.002** (fraction from `IBS.1`; 0.002 from `IBS.2`)

Collapsed to a tree: one source carries <5% of the ancestry, so this graph is behaving as its no-admixture special case.

## Effective sizes shared by IBD and SNP (haploid)

| node | Ne | sd of log Ne |
|---|---|---|
| `EAS` | 1,505,961 | 0.05 |
| `IBS` | 232,399 | 0.06 |
| `TSI` | 1,498,244 | 0.05 |
| `IBS.1` | 231,490 | 0.06 |
| `IBS.2` | 1,497,751 | 0.05 |
| `n1` | 1,501,561 | 0.05 |
| `n2` | 1,498,939 | 0.05 |
| `root` | 11,139 | 0.26 |

log-Ne random-walk step scale tau = 0.808

## Fit quality by component

| component | log-likelihood | n terms | chi2/n |
|---|---|---|---|
| IBD | -7,483.1 | 222 | 99.35 |
| SNP | -39,816.2 | 6 | 13286.60 |

IBD chi2/n uses Palamara model-derived Normal variance; SNP chi2/n uses its block SEs. Values near 1 indicate residuals on the modeled noise scale.

## SNP covariance residuals `(w_hat - W_pred)/w_se`

| | EAS | IBS | TSI |
|---|---|---|---|
| **EAS** | +117.51 | -115.22 | -117.05 |
| **IBS** | -115.22 | +111.70 | +115.48 |
| **TSI** | -117.05 | +115.48 | +114.55 |

![spectrum](spectrum_fit.png)
