# Code and data: *The h-index as an estimator: efficiency, robustness and threshold decisions*

This repository contains the code and the aggregated simulation output for the paper. Running the scripts
regenerates the figures, the tables and the numbers reported in the main text and in the Supplementary
Material. The application's linked candidate-level records contain personal data and are not included;
the scripts of part C rebuild them from public sources.

## Layout

```
R/          scripts; each opens with a Purpose / Produces / Reads header. Run them from the repository root.
data/sim/   aggregated simulation output (read by the report scripts of part B)
data/hill_weissman_sim.rds   cached output of R/hill_weissman_sim.R
figures/    figures of the paper (scripts write here)
```

## Requirements

- R (>= 4.2) with the CRAN packages `ggplot2`, `VGAM`, `statmod` and `jsonlite`:
  `install.packages(c("ggplot2", "VGAM", "statmod", "jsonlite"))`.
- Python 3 with `pandas` (part C only).
- Part B, step 2 runs on a SLURM cluster; the partition names in `R/sim_run.sbatch`, `R/sim_submit.sh` and
  `R/sim_driver.sh` are placeholders to adapt to your site.

## A. Theory simulations (self-contained; no data files, no internet)

| Script | Produces |
|---|---|
| `R/crossing_fig.R` | Figure 1 (`figures/crossing.png`) |
| `R/are_pareto_fig.R` | Supplementary Table S2 and Figure S1 (`figures/are_pareto.png`) |
| `R/are_simulation.R` | Pareto variance-ratio tables and the geometric contrast (console) |
| `R/zeta_corollary_sim.R` | Supplementary Figure S2 (`figures/are_zeta.png`) |
| `R/geometric_corollary_sim.R` | Geometric distribution, Supplementary Remark S2 (console) |
| `R/two_quantity_geometric_check.R` | Geometric plug-in as a function of $n$ and $\sum_i X_i$, Supplementary Proposition S1 (console) |
| `R/winsorized_hill_sim.R` | Supplementary Table S3 and Figures S3 and S4 (`figures/single_outlier_influence.png`, `figures/tradeoff_winsor.png`) |
| `R/hill_weissman_sim.R` | Supplementary Table S4 (writes `data/hill_weissman_sim.rds`) |

## B. Simulation study (Section 6, Supplementary Section S6)

1. `R/sim_calibrate.R`: calibrated parameters of the eight distributions and the scenario list
   (`data/sim/calibration.csv`, Supplementary Table S5; `data/sim/cells.csv`).
2. `R/sim_run.R`, submitted with `bash R/sim_submit.sh <array-range>` (`R/sim_run.sbatch`; `R/sim_driver.sh`
   runs the stopping rule): one task per (scenario, block of 150 fields), written to `data/sim/raw/`
   (not included).
3. `R/sim_aggregate.R`: mean squared error ratios, Monte Carlo standard errors and the stopping rule
   (`data/sim/results.rds`, `data/sim/summary.csv`; a diagnostic plot
   `figures/sim_mse_ratio_diagnostic.png`).
4. `R/sim_decisions.R`: threshold decisions by observed margin (`data/sim/decisions_short_n.csv`).
5. `R/sim_crossing_table.R`: crossover lengths (`data/sim/crossing_table.txt`, Table 2).
6. `R/sim_report.R`: Figure 2 (`figures/sim_mse_ratio.png`), Supplementary Figure S5
   (`figures/sim_mse_ratio_supp.png`), Table 3 and the numbers quoted in Section 6
   (`data/sim/report_numbers.txt`).

Steps 5 and 6 run from the shipped files in `data/sim/`; steps 2 to 4 need the raw output of step 2.

## C. Application: Italian national scientific habilitation (Section 7, Supplementary Section S7)

Download the two public inputs into `data/asn/` (this folder is excluded from the repository):

- `data/asn/github_sosgang/indicatoriASN1618.tsv` from https://github.com/sosgang/coverage_asn
  (`execution/execution_scripts/indicatoriASN1618.tsv`);
- `data/asn/zenodo_4603624/input-anonymized/ND_INPUT_DATA_ANONYMIZED.tsv` from
  https://doi.org/10.5281/zenodo.4603624 (`input-anonymized.zip`).

Then run, in order:

1. `python R/asn_link_deposits.py`: links the two deposits (39,201 applications).
2. `python R/asn_fetch_openalex.py`: retrieves the OpenAlex records of the CV DOIs (set `OPENALEX_API_KEY`
   if you have one; replace the placeholder `MAILTO` address with your own).
3. `python R/asn_build_records.py`: one citation record per application, reconstructed to the deadline.
4. `python R/asn_person_clusters.py`: applications grouped by person (DOI overlap).
5. `Rscript R/asn_analysis.R prorata`, then `lower` and `upper`: estimators for each application.
6. `Rscript R/asn_summarize.R`: Table 4 and Supplementary Table S6.

OpenAlex citation counts grow over time; the paper used responses retrieved in September 2026, so a new
retrieval gives slightly different numbers. The indicator table contains candidate names, and the derived
files carry identifiers that link back to them: keep `data/asn/` local.

## Reproducibility

- The theory scripts use `set.seed(1729)`. The simulation gives field `f` of scenario `i` the L'Ecuyer-CMRG
  stream `(i - 1) * 600 + f`, reached from `set.seed(20260930)`, and stores each field's seed in its output.
  The application's bootstrap uses `set.seed(20260930)`.
- Replicates: B = 8000 for the Pareto, zeta and single-count simulations; B = 4000 for the
  Hill--Weissman simulation; 150 to 600 fields of 30 records per simulation scenario; B = 2000 bootstrap
  resamples in the application.

## License and citation

Code released under the MIT License (see `LICENSE`). For citation details see `CITATION.cff`.
