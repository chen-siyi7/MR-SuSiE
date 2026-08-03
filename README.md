# MRSuSiE

`MRSuSiE` implements joint selection of causal exposures and variants with
direct outcome effects from multivariable Mendelian randomization summary data.

## Installation

```r
install.packages("remotes")
remotes::install_github("chen-siyi7/MR-SuSiE")
```

For a local source archive:

```r
install.packages("MRSuSiE_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Main functions

- `mv_MR_susie_sel_joint()` analyzes continuous-outcome summary statistics.
- `mv_MR_susie_sel_joint_score()` uses signed score statistics for a binary outcome.
- `mvmr_ivw_ld_adjusted()` estimates causal coefficients after selection while
  retaining LD-induced nuisance directions.
- `mv_MR_susie_sel_sequential()` provides the ordering comparator used in the manuscript.

See the function help pages for input definitions and returned diagnostics.
