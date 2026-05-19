# inferMM-enzyme-repro

Reproducible code, input data, and saved output files for the manuscript
_Variance-Aware Estimation and Inference for Michaelis-Menten Regression with Heteroscedastic and Clustered Enzyme-Kinetic Data_.

This repository is a journal-neutral companion repository for the `inferMM`
package and supersedes the earlier `inferMM-cils-repro` naming. The package
repository contains the reusable software interface, while this repository
contains manuscript-oriented simulation scripts, real-data analysis scripts,
raw and processed data, and saved output tables and figures used in the paper.
The manuscript source and PDF files are maintained separately and are not
included here.

## Repository layout

- `scripts/`
  - `michaelis_menten_simulation.R`: main simulation code for the single-curve
    and clustered benchmarking results
  - `bootstrap_ci_sensitivity.R`: small-sample bootstrap sensitivity analyses
  - `real_data.Rmd`: SDL and Alves real-data analyses and figure generation
- `data/raw/`
  - original SDL and Alves input files used in the analyses
- `data/processed/`
  - processed SDL reaction-rate data derived from the raw SDL file
- `results/`
  - `chapter3_simulation_n50/`: main simulation outputs used in the current
    manuscript draft
  - `chapter3_simulation/`: auxiliary simulation outputs, including denser-grid
    and skewed-error summaries
  - `chapter4_real_data_analysis/`: real-data figures and tables

## Dependency

The scripts rely on the companion R package `inferMM`, maintained separately in
its own repository:

```r
# install.packages("remotes")
remotes::install_github("mijeong-kim/inferMM")
```

## Running the analyses

Run commands from the repository root.

- Bootstrap sensitivity analysis:
  - `Rscript scripts/bootstrap_ci_sensitivity.R`
- Real-data analysis and figure generation:
  - `Rscript -e 'rmarkdown::render("scripts/real_data.Rmd")'`
- Main simulation workflow:
  - source `scripts/michaelis_menten_simulation.R` in R and call
    `run_michaelis_menten_simulation(n_obs = 50, use_bayes = FALSE, ...)`
    with an explicit `save_dir`, for example
    `results/chapter3_simulation_n50/`

## Notes

- The current manuscript uses the `n = 50` simulation outputs in
  `results/chapter3_simulation_n50/` as its main benchmark.
- Supplementary analyses include skewed-error robustness summaries,
  design-density checks under the original `n = 100` fixed grid, and
  bootstrap sensitivity analyses at `n = 20` and `n = 30`.
- The real-data analysis script reflects the current package workflow,
  including grouped SDL analyses and pooled Alves comparisons built around
  `group_mm()`.
