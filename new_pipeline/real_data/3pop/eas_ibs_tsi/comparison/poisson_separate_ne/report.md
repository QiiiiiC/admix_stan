# Poisson, separate Ne topology ranking

Populations: **EAS / IBS / TSI**

The weight is a softmax of Pathfinder ELBOs, not a posterior model probability.
Each topology is screened from dispersed MAP starts; distinct high-MAP modes are then fitted independently with Pathfinder. Search diagnostics are retained in `fit.json` and `elbo_table.csv`.
`f` is the fitted ancestry fraction from the first named source branch; tree topologies have no fraction. The closest tree is shown when `min(f, 1-f) < 0.01`.

| rank | topology | f | closest non-admix tree | ELBO | dELBO | ELBO weight |
|---|---|---:|---|---:|---:|---:|
| 1 | [`((EAS,IBS.1),(IBS.2,TSI))`](../../15_((EAS,IBS.1),(IBS.2,TSI))/poisson_separate_ne/report.md) | 0.00548 | `[03] ((IBS,TSI),EAS)` | -499.8 | +0.0 | 1.0000 |
| 2 | [`(((EAS.1,IBS),TSI),EAS.2)`](../../05_(((EAS.1,IBS),TSI),EAS.2)/poisson_separate_ne/report.md) | 0.99851 | `[03] ((IBS,TSI),EAS)` | -529.3 | -29.5 | 0.0000 |
| 3 | [`(((IBS,TSI),EAS.1),EAS.2)`](../../08_(((IBS,TSI),EAS.1),EAS.2)/poisson_separate_ne/report.md) | 0.04199 | `-` | -822.4 | -322.6 | 0.0000 |
| 4 | [`(((EAS.1,TSI),IBS),EAS.2)`](../../07_(((EAS.1,TSI),IBS),EAS.2)/poisson_separate_ne/report.md) | 0.99639 | `[03] ((IBS,TSI),EAS)` | -828.0 | -328.2 | 0.0000 |
| 5 | [`((EAS,TSI.1),(IBS,TSI.2))`](../../21_((EAS,TSI.1),(IBS,TSI.2))/poisson_separate_ne/report.md) | 0.00075 | `[03] ((IBS,TSI),EAS)` | -830.0 | -330.2 | 0.0000 |
| 6 | [`(((EAS,IBS.1),TSI),IBS.2)`](../../11_(((EAS,IBS.1),TSI),IBS.2)/poisson_separate_ne/report.md) | 0.99986 | `[02] ((EAS,TSI),IBS)` | -3687.6 | -3187.8 | 0.0000 |
| 7 | [`(((EAS,IBS.1),IBS.2),TSI)`](../../10_(((EAS,IBS.1),IBS.2),TSI)/poisson_separate_ne/report.md) | 0.99969 | `[01] ((EAS,IBS),TSI)` | -3721.1 | -3221.3 | 0.0000 |
| 8 | [`((EAS.1,IBS),(EAS.2,TSI))`](../../09_((EAS.1,IBS),(EAS.2,TSI))/poisson_separate_ne/report.md) | 0.00000 | `[02] ((EAS,TSI),IBS)` | -3725.6 | -3225.8 | 0.0000 |
| 9 | [`(((EAS.1,IBS),EAS.2),TSI)`](../../04_(((EAS.1,IBS),EAS.2),TSI)/poisson_separate_ne/report.md) | 0.99998 | `[01] ((EAS,IBS),TSI)` | -3733.3 | -3233.5 | 0.0000 |
| 10 | [`(((EAS,TSI.1),IBS),TSI.2)`](../../17_(((EAS,TSI.1),IBS),TSI.2)/poisson_separate_ne/report.md) | 0.99983 | `[01] ((EAS,IBS),TSI)` | -4245.1 | -3745.3 | 0.0000 |
| 11 | [`(((EAS,TSI.1),TSI.2),IBS)`](../../18_(((EAS,TSI.1),TSI.2),IBS)/poisson_separate_ne/report.md) | 0.99948 | `[02] ((EAS,TSI),IBS)` | -4292.3 | -3792.5 | 0.0000 |
| 12 | [`(((EAS.1,TSI),EAS.2),IBS)`](../../06_(((EAS.1,TSI),EAS.2),IBS)/poisson_separate_ne/report.md) | 0.99978 | `[02] ((EAS,TSI),IBS)` | -4298.6 | -3798.8 | 0.0000 |
| 13 | [`(((IBS,TSI.1),EAS),TSI.2)`](../../19_(((IBS,TSI.1),EAS),TSI.2)/poisson_separate_ne/report.md) | 0.81266 | `-` | -4767.8 | -4268.0 | 0.0000 |
| 14 | [`(((IBS.1,TSI),EAS),IBS.2)`](../../13_(((IBS.1,TSI),EAS),IBS.2)/poisson_separate_ne/report.md) | 0.68301 | `-` | -5116.7 | -4616.8 | 0.0000 |
| 15 | [`(((IBS.1,TSI),IBS.2),EAS)`](../../14_(((IBS.1,TSI),IBS.2),EAS)/poisson_separate_ne/report.md) | 0.99834 | `[03] ((IBS,TSI),EAS)` | -5515.2 | -5015.4 | 0.0000 |
| 16 | [`(((IBS,TSI.1),TSI.2),EAS)`](../../20_(((IBS,TSI.1),TSI.2),EAS)/poisson_separate_ne/report.md) | 0.99866 | `[03] ((IBS,TSI),EAS)` | -5526.6 | -5026.8 | 0.0000 |
| 17 | [`((IBS,TSI),EAS)`](../../03_((IBS,TSI),EAS)/poisson_separate_ne/report.md) | - | `-` | -5779.4 | -5279.6 | 0.0000 |
| 18 | [`(((EAS,TSI),IBS.1),IBS.2)`](../../12_(((EAS,TSI),IBS.1),IBS.2)/poisson_separate_ne/report.md) | 0.41937 | `-` | -7156.3 | -6656.5 | 0.0000 |
| 19 | [`(((EAS,IBS),TSI.1),TSI.2)`](../../16_(((EAS,IBS),TSI.1),TSI.2)/poisson_separate_ne/report.md) | 0.44713 | `-` | -7600.6 | -7100.8 | 0.0000 |
| 20 | [`((EAS,IBS),TSI)`](../../01_((EAS,IBS),TSI)/poisson_separate_ne/report.md) | - | `-` | -7618.2 | -7118.4 | 0.0000 |
| 21 | [`((EAS,TSI),IBS)`](../../02_((EAS,TSI),IBS)/poisson_separate_ne/report.md) | - | `-` | -7683.7 | -7183.9 | 0.0000 |

![ranking](elbo_ranking.png)
