# Topology selection on IBS / TSI / CEU

21 models (3 trees + 18 one-admixture graphs), `mixed_model_Nsmooth`, Pathfinder with 6+ seeds, best ELBO kept.

## Ranking

| rank | topology | adm | ELBO | dELBO | restart range | ESS | chi2/n IBD | chi2/n SNP | f | P_model | P_argmax |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [`(((CEU,IBS.1),TSI),IBS.2)`](../05_(((CEU,IBS.1),TSI),IBS.2)/report.md) | 1 | -475.0 | +0.0 | 2.9 | 9.6 | 1.88 | 1.01 | 0.595 | 1.000 | 1.000 |
| 2 | [`(((CEU,TSI.1),IBS),TSI.2)`](../11_(((CEU,TSI.1),IBS),TSI.2)/report.md) | 1 | -566.7 | -91.7 | 0.6 | 4.0 | 2.35 | 10.29 | 0.879 | 0.000 | 0.000 |
| 3 | [`(((CEU,IBS.1),IBS.2),TSI)`](../04_(((CEU,IBS.1),IBS.2),TSI)/report.md) | 1 | -574.9 | -99.9 | 0.9 | 3.3 | 2.87 | 1.93 | 0.983 | 0.000 | 0.000 |
| 4 | [`((CEU,IBS.1),(IBS.2,TSI))`](../09_((CEU,IBS.1),(IBS.2,TSI))/report.md) | 1 | -626.4 | -151.4 | 1.8 | 7.2 | 3.08 | 3.42 | 0.617 | 0.000 | 0.000 |
| 5 | [`(((CEU.1,IBS),CEU.2),TSI)`](../16_(((CEU.1,IBS),CEU.2),TSI)/report.md) | 1 | -653.0 | -178.0 | 1.9 | 2.6 | 3.18 | 6.18 | 1.000 | 0.000 | 0.000 |
| 6 | [`(((CEU.1,IBS),TSI),CEU.2)`](../17_(((CEU.1,IBS),TSI),CEU.2)/report.md) | 1 | -665.0 | -190.0 | 3.2 | 17.9 | 3.11 | 9.54 | 0.175 | 0.000 | 0.000 |
| 7 | [`((CEU.1,IBS),(CEU.2,TSI))`](../21_((CEU.1,IBS),(CEU.2,TSI))/report.md) | 1 | -872.4 | -397.4 | 3.0 | 10.0 | 4.86 | 2.63 | 0.997 | 0.000 | 0.000 |
| 8 | [`((CEU,TSI.1),(IBS,TSI.2))`](../15_((CEU,TSI.1),(IBS,TSI.2))/report.md) | 1 | -1175.0 | -700.0 | 1.4 | 3.1 | 7.38 | 1.09 | 0.989 | 0.000 | 0.000 |
| 9 | [`(((CEU.1,TSI),CEU.2),IBS)`](../18_(((CEU.1,TSI),CEU.2),IBS)/report.md) | 1 | -1467.6 | -992.6 | 1.8 | 3.0 | 11.32 | 10.99 | 0.996 | 0.000 | 0.000 |
| 10 | [`(((CEU.1,TSI),IBS),CEU.2)`](../19_(((CEU.1,TSI),IBS),CEU.2)/report.md) | 1 | -1478.8 | -1003.8 | 2.2 | 7.4 | 11.40 | 7.28 | 1.000 | 0.000 | 0.000 |
| 11 | [`(((IBS,TSI.1),CEU),TSI.2)`](../13_(((IBS,TSI.1),CEU),TSI.2)/report.md) | 1 | -1535.1 | -1060.1 | 1.4 | 7.0 | 7.56 | 128.35 | 0.446 | 0.000 | 0.000 |
| 12 | [`(((CEU,TSI.1),TSI.2),IBS)`](../12_(((CEU,TSI.1),TSI.2),IBS)/report.md) | 1 | -1897.6 | -1422.6 | 1.8 | 11.1 | 9.63 | 203.33 | 0.999 | 0.000 | 0.000 |
| 13 | [`((CEU,IBS),TSI)`](../01_((CEU,IBS),TSI)/report.md) | 0 | -1920.1 | -1445.1 | 2.2 | 28.9 | 4.95 | 359.47 | - | 0.000 | 0.000 |
| 14 | [`(((CEU,IBS),TSI.1),TSI.2)`](../10_(((CEU,IBS),TSI.1),TSI.2)/report.md) | 1 | -1928.9 | -1453.9 | 5.2 | 6.5 | 4.93 | 347.91 | 0.542 | 0.000 | 0.000 |
| 15 | [`(((IBS,TSI),CEU.1),CEU.2)`](../20_(((IBS,TSI),CEU.1),CEU.2)/report.md) | 1 | -1980.6 | -1505.6 | 2.4 | 11.6 | 8.34 | 249.47 | 0.015 | 0.000 | 0.000 |
| 16 | [`(((IBS.1,TSI),CEU),IBS.2)`](../07_(((IBS.1,TSI),CEU),IBS.2)/report.md) | 1 | -2140.6 | -1665.6 | 2.2 | 1.8 | 25.75 | 11.68 | 0.443 | 0.000 | 0.000 |
| 17 | [`(((IBS.1,TSI),IBS.2),CEU)`](../08_(((IBS.1,TSI),IBS.2),CEU)/report.md) | 1 | -2151.8 | -1676.8 | 2.2 | 2.7 | 27.18 | 13.75 | 0.881 | 0.000 | 0.000 |
| 18 | [`(((IBS,TSI.1),TSI.2),CEU)`](../14_(((IBS,TSI.1),TSI.2),CEU)/report.md) | 1 | -2592.5 | -2117.5 | 3.3 | 8.7 | 22.07 | 184.57 | 0.998 | 0.000 | 0.000 |
| 19 | [`(((CEU,TSI),IBS.1),IBS.2)`](../06_(((CEU,TSI),IBS.1),IBS.2)/report.md) | 1 | -2698.4 | -2223.5 | 5.1 | 7.9 | 10.76 | 399.72 | 0.156 | 0.000 | 0.000 |
| 20 | [`((IBS,TSI),CEU)`](../03_((IBS,TSI),CEU)/report.md) | 0 | -2820.5 | -2345.5 | 2.2 | 1.1 | 27.82 | 208.96 | - | 0.000 | 0.000 |
| 21 | [`((CEU,TSI),IBS)`](../02_((CEU,TSI),IBS)/report.md) | 0 | -2916.3 | -2441.3 | 0.9 | 15.4 | 14.90 | 388.68 | - | 0.000 | 0.000 |

![ranking](elbo_ranking.png)

## Selection

- **Best model: `(((CEU,IBS.1),TSI),IBS.2)`**, ELBO -475.0, P_model = 1.000, P_argmax = 1.000.
- Runner-up `(((CEU,TSI.1),IBS),TSI.2)` is -91.7 nats behind. The winner's 6 restarts span 2.9 nats, and **6 of them beat the runner-up's best restart** -- the lead does not depend on the restart budget.
- Best tree `((CEU,IBS),TSI)` vs best admixture graph `(((CEU,IBS.1),TSI),IBS.2)`: +1445.1 nats for 3 extra Ne, 2 extra times and 1 fraction.

## How to read the two probabilities

`P_model` is softmax over the 21 ELBOs with a uniform model prior -- the posterior model probability. `P_argmax` bootstraps the Pathfinder draws and counts how often each model wins, so it measures only whether the RANKING is stable against Monte Carlo error. Both saturate at 1.000/0.000 whenever the gaps exceed a few nats, which they do here; the informative number is dELBO.

## Caveats that travel with these numbers

1. The ELBO is a lower bound on log Z and the gap KL(q||p) differs per model, so a model whose posterior happens to be more Gaussian is rewarded for that alone. The importance-sampling `logZ` would correct it, but the ESS column shows the weights are degenerate (max ESS 29 of 4000), so `logZ` cannot be trusted here either.
2. An admixture graph splits the admixed leaf into two source branches with separate Ne -- it buys a piecewise Ne, which these data independently want. A win by an admixture graph is evidence for non-constant Ne at least as much as for admixture, and three leaves cannot separate the two. Check each folder's residual row for a monotone tilt before reading any result as admixture.
3. These probabilities are conditional on the 21 listed models being the whole hypothesis space. Two admixture events, and non-constant Ne without admixture, are not in it.
