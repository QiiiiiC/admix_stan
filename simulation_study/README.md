# simulation_study

One folder per simulation study, named after its **true (generating)
topology**.  Every study imports the shared machinery from `../methods/`
and keeps only what is specific to it: the generating history, the candidate
graphs, how candidate nodes map onto true branches, and its own tables and
figures.

Each study folder follows the same layout:

| File | Role |
|---|---|
| `config.json` | every number of the design (times, sizes, samples, blocks, bins, fit settings) |
| `study.py` | generating history, candidates, truth mappings; binds config into `methods` |
| `run_study.py` | resumable stages `plan / compile / simulate / fit / summarize / calibrate`, sliceable across workers |
| `execute_full.py` | supervisor: compile once, simulate once, fit in parallel, summarise |
| `visualize.py`, `plot_*.py` | tables, report and figures (study-specific) |
| `test_design.py` | design invariants |
| `README.md`, `VALIDATION.md` | what the study asks and how it was checked |
| `summary/` | tracked results: CSVs, `report.md`, `figures/` |
| `runs/` | simulated blocks and fits (git-ignored, large) |

## Studies

| Folder | True topology | Question | Status |
|---|---|---|---|
| `growth_asymmetric_loop` | 3 leaves: recent loop in b (10–25 gen) inside a b admixture (30 gen); growth on a (10×) and b (4×), 3:1 size asymmetries | Does a separate IBD/SNP Ne help when every candidate is misspecified in Ne? | Complete (ported from `new_pipeline/Nfixed/`, verified equivalent) |
| `admix_thru_b_plus_ancient` | 4 leaves, Ne 10,000: recent admix-through-b loop on a (20→60→100 gen, α = 0.5) + ancient 70/30 admixture on c (500 gen) | Does the mixed model beat IBD-only and SNP-only on both events, and from what genome length (50–1000 cM)? | Complete (12,000 fits, 0 failed; rerun of `new_pipeline/Nfixed/topology_admix_thru_b_plus_ancient` with Poisson IBD) |
| `four_leaf_b_admixture_varying_ne` | 4 leaves: b admixed at 20 gen (0.7 from the a side), merges 60/100/200/400; a different Ne on every branch (5,000–20,000) | Parameter recovery only: which of the 16 parameters do IBD-only, Mixed (Nvarying) and Mixed (Nsmooth) recover, and from what length? | Complete (4,800 fits, 0 failed); rerun through the growth pipeline in `runs/pipeline` with true IBD and hap-IBD: Pathfinder complete (6,400 fits, 0 failed), NUTS running |
| `four_leaf_b_admixture_growth_ne` | Same graph; growth-like Ne (recent branches 23k–122k, root 10k); 22 bins from 1 cM; true IBD and hap-IBD | Mixed vs IBD-only Nsmooth: T_true vs tree (identification), parameter recovery, spectrum fit at 1000 cM, exact log Z vs best-path ELBO | Pathfinder complete (6,400 fits, 0 failed); NUTS (160, true IBD) running |
