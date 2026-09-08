# Normal, shared Ne, grid topology ranking

Populations: **EAS / IBS / TSI**

The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.
Each topology is screened from dispersed MAP starts; distinct high-MAP modes are then fitted independently with Pathfinder. Search diagnostics are retained in `fit.json` and `elbo_table.csv`.
`f` is the fitted ancestry fraction from the first named source branch; tree topologies have no fraction. The closest tree is shown when `min(f, 1-f) < 0.01`.

| rank | topology | f | closest non-admix tree | ELBO | dELBO | ELBO weight |
|---|---|---:|---|---:|---:|---:|
| 1 | [`(((EAS.1,IBS),TSI),EAS.2)`](../../05_(((EAS.1,IBS),TSI),EAS.2)/normal_grid_shared_ne/report.md) | 0.99222 | `[03] ((IBS,TSI),EAS)` | +3994.5 | +0.0 | 1.0000 |
| 2 | [`(((IBS,TSI),EAS.1),EAS.2)`](../../08_(((IBS,TSI),EAS.1),EAS.2)/normal_grid_shared_ne/report.md) | 0.80572 | `-` | +3906.6 | -88.0 | 0.0000 |
| 3 | [`(((EAS.1,TSI),IBS),EAS.2)`](../../07_(((EAS.1,TSI),IBS),EAS.2)/normal_grid_shared_ne/report.md) | 0.98690 | `-` | +3813.6 | -180.9 | 0.0000 |
| 4 | [`((EAS,IBS.1),(IBS.2,TSI))`](../../15_((EAS,IBS.1),(IBS.2,TSI))/normal_grid_shared_ne/report.md) | 0.00379 | `[03] ((IBS,TSI),EAS)` | +3098.2 | -896.3 | 0.0000 |
| 5 | [`((EAS,TSI.1),(IBS,TSI.2))`](../../21_((EAS,TSI.1),(IBS,TSI.2))/normal_grid_shared_ne/report.md) | 0.00100 | `[03] ((IBS,TSI),EAS)` | +2890.8 | -1103.7 | 0.0000 |
| 6 | [`(((IBS.1,TSI),IBS.2),EAS)`](../../14_(((IBS.1,TSI),IBS.2),EAS)/normal_grid_shared_ne/report.md) | 0.11719 | `-` | +383.9 | -3610.6 | 0.0000 |
| 7 | [`(((IBS,TSI.1),TSI.2),EAS)`](../../20_(((IBS,TSI.1),TSI.2),EAS)/normal_grid_shared_ne/report.md) | 0.85723 | `-` | +285.8 | -3708.8 | 0.0000 |
| 8 | [`(((IBS.1,TSI),EAS),IBS.2)`](../../13_(((IBS.1,TSI),EAS),IBS.2)/normal_grid_shared_ne/report.md) | 0.00227 | `[03] ((IBS,TSI),EAS)` | -1093.9 | -5088.4 | 0.0000 |
| 9 | [`((IBS,TSI),EAS)`](../../03_((IBS,TSI),EAS)/normal_grid_shared_ne/report.md) | - | `-` | -1104.6 | -5099.1 | 0.0000 |
| 10 | [`(((IBS,TSI.1),EAS),TSI.2)`](../../19_(((IBS,TSI.1),EAS),TSI.2)/normal_grid_shared_ne/report.md) | 0.00006 | `[03] ((IBS,TSI),EAS)` | -1120.7 | -5115.3 | 0.0000 |
| 11 | [`((EAS.1,IBS),(EAS.2,TSI))`](../../09_((EAS.1,IBS),(EAS.2,TSI))/normal_grid_shared_ne/report.md) | 0.99950 | `[01] ((EAS,IBS),TSI)` | -1455.8 | -5450.3 | 0.0000 |
| 12 | [`(((EAS,IBS.1),IBS.2),TSI)`](../../10_(((EAS,IBS.1),IBS.2),TSI)/normal_grid_shared_ne/report.md) | 0.99906 | `[01] ((EAS,IBS),TSI)` | -35317.3 | -39311.9 | 0.0000 |
| 13 | [`(((EAS,IBS.1),TSI),IBS.2)`](../../11_(((EAS,IBS.1),TSI),IBS.2)/normal_grid_shared_ne/report.md) | 0.99979 | `[02] ((EAS,TSI),IBS)` | -35320.8 | -39315.3 | 0.0000 |
| 14 | [`(((EAS.1,IBS),EAS.2),TSI)`](../../04_(((EAS.1,IBS),EAS.2),TSI)/normal_grid_shared_ne/report.md) | 0.99899 | `[01] ((EAS,IBS),TSI)` | -35359.8 | -39354.3 | 0.0000 |
| 15 | [`(((EAS,TSI.1),IBS),TSI.2)`](../../17_(((EAS,TSI.1),IBS),TSI.2)/normal_grid_shared_ne/report.md) | 0.99901 | `[01] ((EAS,IBS),TSI)` | -35481.3 | -39475.9 | 0.0000 |
| 16 | [`(((EAS,TSI.1),TSI.2),IBS)`](../../18_(((EAS,TSI.1),TSI.2),IBS)/normal_grid_shared_ne/report.md) | 0.99825 | `[02] ((EAS,TSI),IBS)` | -35538.0 | -39532.6 | 0.0000 |
| 17 | [`(((EAS.1,TSI),EAS.2),IBS)`](../../06_(((EAS.1,TSI),EAS.2),IBS)/normal_grid_shared_ne/report.md) | 0.99988 | `[02] ((EAS,TSI),IBS)` | -36501.9 | -40496.5 | 0.0000 |
| 18 | [`(((EAS,TSI),IBS.1),IBS.2)`](../../12_(((EAS,TSI),IBS.1),IBS.2)/normal_grid_shared_ne/report.md) | 0.03479 | `-` | -37200.1 | -41194.6 | 0.0000 |
| 19 | [`(((EAS,IBS),TSI.1),TSI.2)`](../../16_(((EAS,IBS),TSI.1),TSI.2)/normal_grid_shared_ne/report.md) | 0.02942 | `-` | -37391.5 | -41386.0 | 0.0000 |
| 20 | [`((EAS,IBS),TSI)`](../../01_((EAS,IBS),TSI)/normal_grid_shared_ne/report.md) | - | `-` | -40954.1 | -44948.6 | 0.0000 |
| 21 | [`((EAS,TSI),IBS)`](../../02_((EAS,TSI),IBS)/normal_grid_shared_ne/report.md) | - | `-` | -40978.0 | -44972.5 | 0.0000 |

![ranking](elbo_ranking.png)
