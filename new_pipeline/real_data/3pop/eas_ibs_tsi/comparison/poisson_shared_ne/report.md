# Poisson, shared Ne topology ranking

Populations: **EAS / IBS / TSI**

The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.

| rank | topology | ELBO | dELBO | ELBO weight |
|---|---|---:|---:|---:|
| 1 | [`((EAS,IBS.1),(IBS.2,TSI))`](../../15_((EAS,IBS.1),(IBS.2,TSI))/poisson_shared_ne/report.md) | -1013.0 | +0.0 | 1.0000 |
| 2 | [`((EAS,TSI.1),(IBS,TSI.2))`](../../21_((EAS,TSI.1),(IBS,TSI.2))/poisson_shared_ne/report.md) | -1058.3 | -45.3 | 0.0000 |
| 3 | [`(((EAS,TSI.1),IBS),TSI.2)`](../../17_(((EAS,TSI.1),IBS),TSI.2)/poisson_shared_ne/report.md) | -5114.2 | -4101.2 | 0.0000 |
| 4 | [`(((EAS,IBS.1),IBS.2),TSI)`](../../10_(((EAS,IBS.1),IBS.2),TSI)/poisson_shared_ne/report.md) | -5132.6 | -4119.6 | 0.0000 |
| 5 | [`(((EAS,TSI.1),TSI.2),IBS)`](../../18_(((EAS,TSI.1),TSI.2),IBS)/poisson_shared_ne/report.md) | -5318.3 | -4305.3 | 0.0000 |
| 6 | [`(((EAS,IBS.1),TSI),IBS.2)`](../../11_(((EAS,IBS.1),TSI),IBS.2)/poisson_shared_ne/report.md) | -5633.3 | -4620.3 | 0.0000 |
| 7 | [`((IBS,TSI),EAS)`](../../03_((IBS,TSI),EAS)/poisson_shared_ne/report.md) | -5922.3 | -4909.3 | 0.0000 |
| 8 | [`(((IBS,TSI.1),EAS),TSI.2)`](../../19_(((IBS,TSI.1),EAS),TSI.2)/poisson_shared_ne/report.md) | -5979.1 | -4966.1 | 0.0000 |
| 9 | [`(((EAS.1,TSI),IBS),EAS.2)`](../../07_(((EAS.1,TSI),IBS),EAS.2)/poisson_shared_ne/report.md) | -6052.5 | -5039.4 | 0.0000 |
| 10 | [`(((IBS,TSI),EAS.1),EAS.2)`](../../08_(((IBS,TSI),EAS.1),EAS.2)/poisson_shared_ne/report.md) | -6057.3 | -5044.3 | 0.0000 |
| 11 | [`(((IBS,TSI.1),TSI.2),EAS)`](../../20_(((IBS,TSI.1),TSI.2),EAS)/poisson_shared_ne/report.md) | -6113.6 | -5100.6 | 0.0000 |
| 12 | [`(((IBS.1,TSI),EAS),IBS.2)`](../../13_(((IBS.1,TSI),EAS),IBS.2)/poisson_shared_ne/report.md) | -6115.0 | -5101.9 | 0.0000 |
| 13 | [`(((EAS.1,IBS),TSI),EAS.2)`](../../05_(((EAS.1,IBS),TSI),EAS.2)/poisson_shared_ne/report.md) | -6163.6 | -5150.6 | 0.0000 |
| 14 | [`(((IBS.1,TSI),IBS.2),EAS)`](../../14_(((IBS.1,TSI),IBS.2),EAS)/poisson_shared_ne/report.md) | -6212.0 | -5199.0 | 0.0000 |
| 15 | [`(((EAS.1,IBS),EAS.2),TSI)`](../../04_(((EAS.1,IBS),EAS.2),TSI)/poisson_shared_ne/report.md) | -42695.7 | -41682.6 | 0.0000 |
| 16 | [`(((EAS.1,TSI),EAS.2),IBS)`](../../06_(((EAS.1,TSI),EAS.2),IBS)/poisson_shared_ne/report.md) | -42821.6 | -41808.6 | 0.0000 |
| 17 | [`((EAS.1,IBS),(EAS.2,TSI))`](../../09_((EAS.1,IBS),(EAS.2,TSI))/poisson_shared_ne/report.md) | -43306.7 | -42293.7 | 0.0000 |
| 18 | [`(((EAS,TSI),IBS.1),IBS.2)`](../../12_(((EAS,TSI),IBS.1),IBS.2)/poisson_shared_ne/report.md) | -43490.9 | -42477.9 | 0.0000 |
| 19 | [`((EAS,IBS),TSI)`](../../01_((EAS,IBS),TSI)/poisson_shared_ne/report.md) | -43568.2 | -42555.2 | 0.0000 |
| 20 | [`(((EAS,IBS),TSI.1),TSI.2)`](../../16_(((EAS,IBS),TSI.1),TSI.2)/poisson_shared_ne/report.md) | -43686.0 | -42673.0 | 0.0000 |
| 21 | [`((EAS,TSI),IBS)`](../../02_((EAS,TSI),IBS)/poisson_shared_ne/report.md) | -43992.1 | -42979.1 | 0.0000 |

![ranking](elbo_ranking.png)
