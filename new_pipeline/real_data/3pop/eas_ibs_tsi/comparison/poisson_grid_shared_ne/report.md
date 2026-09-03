# Poisson, shared Ne, grid topology ranking

Populations: **EAS / IBS / TSI**

The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.

| rank | topology | ELBO | dELBO | ELBO weight |
|---|---|---:|---:|---:|
| 1 | [`((EAS,IBS.1),(IBS.2,TSI))`](../../15_((EAS,IBS.1),(IBS.2,TSI))/poisson_grid_shared_ne/report.md) | -794.7 | +0.0 | 1.0000 |
| 2 | [`((EAS,TSI.1),(IBS,TSI.2))`](../../21_((EAS,TSI.1),(IBS,TSI.2))/poisson_grid_shared_ne/report.md) | -889.1 | -94.4 | 0.0000 |
| 3 | [`(((EAS,TSI.1),IBS),TSI.2)`](../../17_(((EAS,TSI.1),IBS),TSI.2)/poisson_grid_shared_ne/report.md) | -4821.0 | -4026.3 | 0.0000 |
| 4 | [`(((EAS,IBS.1),IBS.2),TSI)`](../../10_(((EAS,IBS.1),IBS.2),TSI)/poisson_grid_shared_ne/report.md) | -4823.3 | -4028.5 | 0.0000 |
| 5 | [`(((EAS,IBS.1),TSI),IBS.2)`](../../11_(((EAS,IBS.1),TSI),IBS.2)/poisson_grid_shared_ne/report.md) | -5019.7 | -4225.0 | 0.0000 |
| 6 | [`(((IBS,TSI),EAS.1),EAS.2)`](../../08_(((IBS,TSI),EAS.1),EAS.2)/poisson_grid_shared_ne/report.md) | -5161.2 | -4366.4 | 0.0000 |
| 7 | [`(((EAS,TSI.1),TSI.2),IBS)`](../../18_(((EAS,TSI.1),TSI.2),IBS)/poisson_grid_shared_ne/report.md) | -5181.2 | -4386.4 | 0.0000 |
| 8 | [`(((EAS.1,IBS),TSI),EAS.2)`](../../05_(((EAS.1,IBS),TSI),EAS.2)/poisson_grid_shared_ne/report.md) | -5210.5 | -4415.8 | 0.0000 |
| 9 | [`((IBS,TSI),EAS)`](../../03_((IBS,TSI),EAS)/poisson_grid_shared_ne/report.md) | -5249.7 | -4455.0 | 0.0000 |
| 10 | [`(((EAS.1,TSI),IBS),EAS.2)`](../../07_(((EAS.1,TSI),IBS),EAS.2)/poisson_grid_shared_ne/report.md) | -5272.7 | -4478.0 | 0.0000 |
| 11 | [`(((IBS.1,TSI),EAS),IBS.2)`](../../13_(((IBS.1,TSI),EAS),IBS.2)/poisson_grid_shared_ne/report.md) | -5302.1 | -4507.4 | 0.0000 |
| 12 | [`(((IBS,TSI.1),TSI.2),EAS)`](../../20_(((IBS,TSI.1),TSI.2),EAS)/poisson_grid_shared_ne/report.md) | -5347.0 | -4552.3 | 0.0000 |
| 13 | [`(((IBS.1,TSI),IBS.2),EAS)`](../../14_(((IBS.1,TSI),IBS.2),EAS)/poisson_grid_shared_ne/report.md) | -5353.4 | -4558.6 | 0.0000 |
| 14 | [`(((IBS,TSI.1),EAS),TSI.2)`](../../19_(((IBS,TSI.1),EAS),TSI.2)/poisson_grid_shared_ne/report.md) | -5417.3 | -4622.5 | 0.0000 |
| 15 | [`(((EAS.1,IBS),EAS.2),TSI)`](../../04_(((EAS.1,IBS),EAS.2),TSI)/poisson_grid_shared_ne/report.md) | -7421.8 | -6627.1 | 0.0000 |
| 16 | [`((EAS.1,IBS),(EAS.2,TSI))`](../../09_((EAS.1,IBS),(EAS.2,TSI))/poisson_grid_shared_ne/report.md) | -42151.6 | -41356.8 | 0.0000 |
| 17 | [`(((EAS.1,TSI),EAS.2),IBS)`](../../06_(((EAS.1,TSI),EAS.2),IBS)/poisson_grid_shared_ne/report.md) | -42160.8 | -41366.1 | 0.0000 |
| 18 | [`(((EAS,TSI),IBS.1),IBS.2)`](../../12_(((EAS,TSI),IBS.1),IBS.2)/poisson_grid_shared_ne/report.md) | -43042.7 | -42248.0 | 0.0000 |
| 19 | [`((EAS,IBS),TSI)`](../../01_((EAS,IBS),TSI)/poisson_grid_shared_ne/report.md) | -43178.3 | -42383.6 | 0.0000 |
| 20 | [`(((EAS,IBS),TSI.1),TSI.2)`](../../16_(((EAS,IBS),TSI.1),TSI.2)/poisson_grid_shared_ne/report.md) | -43238.5 | -42443.8 | 0.0000 |
| 21 | [`((EAS,TSI),IBS)`](../../02_((EAS,TSI),IBS)/poisson_grid_shared_ne/report.md) | -43444.4 | -42649.7 | 0.0000 |

![ranking](elbo_ranking.png)
