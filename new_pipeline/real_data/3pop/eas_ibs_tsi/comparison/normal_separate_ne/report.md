# Normal, separate Ne topology ranking

Populations: **EAS / IBS / TSI**

The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.

| rank | topology | ELBO | dELBO | ELBO weight |
|---|---|---:|---:|---:|
| 1 | [`(((IBS,TSI),EAS.1),EAS.2)`](../../08_(((IBS,TSI),EAS.1),EAS.2)/normal_separate_ne/report.md) | +3815.4 | +0.0 | 1.0000 |
| 2 | [`(((EAS.1,IBS),TSI),EAS.2)`](../../05_(((EAS.1,IBS),TSI),EAS.2)/normal_separate_ne/report.md) | +696.1 | -3119.3 | 0.0000 |
| 3 | [`(((IBS,TSI.1),EAS),TSI.2)`](../../19_(((IBS,TSI.1),EAS),TSI.2)/normal_separate_ne/report.md) | +646.9 | -3168.5 | 0.0000 |
| 4 | [`(((EAS,IBS.1),IBS.2),TSI)`](../../10_(((EAS,IBS.1),IBS.2),TSI)/normal_separate_ne/report.md) | +436.2 | -3379.3 | 0.0000 |
| 5 | [`(((EAS.1,IBS),EAS.2),TSI)`](../../04_(((EAS.1,IBS),EAS.2),TSI)/normal_separate_ne/report.md) | +302.6 | -3512.8 | 0.0000 |
| 6 | [`((EAS.1,IBS),(EAS.2,TSI))`](../../09_((EAS.1,IBS),(EAS.2,TSI))/normal_separate_ne/report.md) | -23.1 | -3838.6 | 0.0000 |
| 7 | [`(((EAS.1,TSI),EAS.2),IBS)`](../../06_(((EAS.1,TSI),EAS.2),IBS)/normal_separate_ne/report.md) | -27.3 | -3842.8 | 0.0000 |
| 8 | [`(((IBS.1,TSI),IBS.2),EAS)`](../../14_(((IBS.1,TSI),IBS.2),EAS)/normal_separate_ne/report.md) | -73.2 | -3888.6 | 0.0000 |
| 9 | [`(((IBS.1,TSI),EAS),IBS.2)`](../../13_(((IBS.1,TSI),EAS),IBS.2)/normal_separate_ne/report.md) | -83.9 | -3899.4 | 0.0000 |
| 10 | [`(((IBS,TSI.1),TSI.2),EAS)`](../../20_(((IBS,TSI.1),TSI.2),EAS)/normal_separate_ne/report.md) | -97.1 | -3912.6 | 0.0000 |
| 11 | [`((IBS,TSI),EAS)`](../../03_((IBS,TSI),EAS)/normal_separate_ne/report.md) | -334.3 | -4149.8 | 0.0000 |
| 12 | [`(((EAS,TSI.1),IBS),TSI.2)`](../../17_(((EAS,TSI.1),IBS),TSI.2)/normal_separate_ne/report.md) | -493.1 | -4308.6 | 0.0000 |
| 13 | [`(((EAS,TSI.1),TSI.2),IBS)`](../../18_(((EAS,TSI.1),TSI.2),IBS)/normal_separate_ne/report.md) | -709.4 | -4524.8 | 0.0000 |
| 14 | [`((EAS,IBS),TSI)`](../../01_((EAS,IBS),TSI)/normal_separate_ne/report.md) | -1566.0 | -5381.5 | 0.0000 |
| 15 | [`(((EAS,IBS),TSI.1),TSI.2)`](../../16_(((EAS,IBS),TSI.1),TSI.2)/normal_separate_ne/report.md) | -1595.4 | -5410.8 | 0.0000 |
| 16 | [`(((EAS.1,TSI),IBS),EAS.2)`](../../07_(((EAS.1,TSI),IBS),EAS.2)/normal_separate_ne/report.md) | -3172.1 | -6987.5 | 0.0000 |
| 17 | [`((EAS,IBS.1),(IBS.2,TSI))`](../../15_((EAS,IBS.1),(IBS.2,TSI))/normal_separate_ne/report.md) | -3191.8 | -7007.2 | 0.0000 |
| 18 | [`((EAS,TSI.1),(IBS,TSI.2))`](../../21_((EAS,TSI.1),(IBS,TSI.2))/normal_separate_ne/report.md) | -3197.8 | -7013.2 | 0.0000 |
| 19 | [`(((EAS,IBS.1),TSI),IBS.2)`](../../11_(((EAS,IBS.1),TSI),IBS.2)/normal_separate_ne/report.md) | -26035.2 | -29850.6 | 0.0000 |
| 20 | [`(((EAS,TSI),IBS.1),IBS.2)`](../../12_(((EAS,TSI),IBS.1),IBS.2)/normal_separate_ne/report.md) | -27838.1 | -31653.6 | 0.0000 |
| 21 | [`((EAS,TSI),IBS)`](../../02_((EAS,TSI),IBS)/normal_separate_ne/report.md) | -28149.3 | -31964.7 | 0.0000 |

![ranking](elbo_ranking.png)
