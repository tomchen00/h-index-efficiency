# Code and data: *Robust but inefficient: the statistical efficiency of threshold bibliometric indices*

This repository contains the code and cached data for the paper. Running the scripts regenerates the
figures, tables, and the numerical summaries reported in it. It holds ten R scripts (each with a
Purpose / Produces / Reads / Requires header), the cached OpenAlex citation-count and publication-year
data for the six-author illustration, and an output folder for the figures.

## Requirements

- R (>= 4.2). The scripts use only CRAN packages.
- Install the packages once:

  ```r
  install.packages(c("ggplot2", "tidyr", "fitdistrplus", "actuar", "VGAM", "jsonlite"))
  ```

Only `openalex_panel.R` needs internet access, and only if you delete the cached data (see below).

## How to run

Set the working directory to this folder (the one containing this README), then run a script with
`Rscript <name>.R` from a shell, or `source("<name>.R")` inside R. Figures are written to `figures/`;
data caches live in `data/`. Each script re-creates those folders if your unzip tool dropped them.

**A. Self-contained theory simulations** (no data files, no internet; run in any order):

| Script | Produces |
|---|---|
| `are_pareto_fig.R` | Table 1 and `figures/are_pareto.png` (simulation, finite-n prediction, and leading rate) |
| `are_simulation.R` | Pareto variance-ratio tables (alpha = 1 and 2) and the geometric (light-tail) contrast (console) |
| `zeta_corollary_sim.R` | Discrete zeta corollary and `figures/are_zeta.png` (curves drawn to n = 1e6) |
| `winsorized_hill_sim.R` | Efficiency-robustness table and `figures/tradeoff_winsor.png`, `figures/single_outlier_influence.png` |
| `geometric_corollary_sim.R` | Light-tailed (geometric) lattice-degeneracy check (console) |
| `two_quantity_geometric_check.R` | Sufficiency check: geometric `h` is a function of `(n, C)` (console) |

**B. OpenAlex illustration** (run in order; later steps read earlier caches):

1. `openalex_panel.R`: fetches the six authors and writes `data/authors_data.rds`, `data/panel_summary.rds`.
   The caches are shipped, so this step is optional; it re-pulls only if `data/authors_data.rds` is missing.
2. `fit_families_full.R`: fits the tail families, writes `data/family_fits.rds`.
3. `fit_families_boot.R`: Burr bootstrap, writes `data/burr_boot.rds` (B = 4000; the most time-consuming step).
4. `openalex_panel_figs.R`: writes `figures/openalex_loglog.png`, `figures/openalex_biasvar.png`.

Because the four `.rds` caches are included, you can run step 4 directly to regenerate the OpenAlex-based
figures without internet or the refits.

## Reproducibility

- Every stochastic script sets `set.seed(1729)`. Replicate counts: simulation figures use B = 8000;
  the Burr bootstrap uses B = 4000.
- `data/authors_data.rds` is the raw OpenAlex pull (`cited_by_count`, `publication_year`) for the six
  authors in the paper, retrieved from the public API (`api.openalex.org`) on 6 June 2026. Re-pulling later
  gives different counts as citations accrue; the cache preserves the OpenAlex snapshot used in the paper.
- `data/{panel_summary, family_fits, burr_boot}.rds` are derived caches; delete them to recompute from
  `authors_data.rds`.
- `openalex_panel.R` sends a placeholder contact email (`anonymous@example.org`) to OpenAlex's polite pool;
  replace it with your own address if you re-run the fetch. It is not needed for the cached run.

## Files

```
README.md
*.R                      nine scripts (each opens with a Purpose / Produces / Reads / Requires header)
data/  authors_data.rds  raw OpenAlex pull (six authors)
       panel_summary.rds  per-author empirical/plug-in summary and bootstrap variance ratios
       family_fits.rds    KS distance and model-implied h per tail family
       burr_boot.rds      Burr bootstrap calibration gap and variance ratio
figures/                 output folder (figures are written here)
```

Data source: OpenAlex (Priem, Piwowar & Orr, 2022), public API, retrieved 6 June 2026.

## License and citation

Code released under the MIT License (see `LICENSE`). For citation details see `CITATION.cff`.
