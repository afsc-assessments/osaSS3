# osaSS3 (experimental)

An installable R package for iterative Stock Synthesis 3 (SS3) multinomial composition reweighting using standard deviations of one-step-ahead (OSA) residuals, with final per-fleet/stratum **four-panel composite PDF diagnostics** based on `afscOSA`.

> **Research software, not validated for production assessments.** The inverse-SDNR-squared update is a heuristic, not a formally validated likelihood-weight estimator. Inspect the residuals, fitted outputs, sample sizes, and model behavior before interpreting the results.

## Install

```r
install.packages("remotes")
# After this package is pushed to a GitHub repository:
# remotes::install_github("OWNER/osaSS3")
```

Install `r4ss`, `afscOSA`, and their dependencies according to their documentation. The SS3 executable is **not** included.

## Use

```r
library(osaSS3)
fit <- osa_reweight_ss3(
  model_dir = "C:/SS3/MyModel",
  plot_final = TRUE,
  max_iter = 8
)
fit$history
fit$final_plot_files
```

The model folder must contain `starter.ss`, the indicated data and control files, forecast input if needed, and a runnable SS3 executable (`ss3.exe` in the folder or a path provided as `ss_exe`).

`OSA_reweight_runs/iter_XX/final_osa_plots/` contains the final diagnostic PDFs. Each stratum gets a four-panel page: aggregated proportions, OSA Q-Q, OSA residual bubbles, and Pearson residual bubbles.

## Scope and limitations

- Detects positive variance-adjustment factors **4 (length)** and **5 (age)** in the control file.
- Fits separate OSA residual analyses by available strata, including age-error groups, then pools residuals within each fleet/factor for a single weight update.
- Does **not** support conditional age-at-length compositions, generalized size data, Dirichlet-multinomial likelihoods, or all model configurations.
- Refuses to overwrite an existing output directory.
- The last *proposed* next weight is not fitted if the loop terminates at the iteration limit.
- This implementation has not been tested end-to-end in R in the authoring environment; report bugs with minimal reproducible examples.

## Development

The package is based on the generic OSA script developed for SS3 Pacific cod composition tuning. The license file uses MIT; replace the placeholder maintainer information in `DESCRIPTION` before public release.
