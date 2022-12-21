
<!-- README.md is generated from README.Rmd. Please edit that file -->

# agriutilities <a href="https://apariciojohan.github.io/agriutilities/"><img src="man/figures/logo.png" align="right" height="132" /></a>

<!-- badges: start -->
<!-- badges: end -->

The goal of this package is to provide some additional functions for the
analysis of field trials in agriculture. It extends the functionalities
of packages such as as SpATS, lme4, ASReml-R, and statgenSTA.

## Installation

You can install the development version of agriutilities from
[GitHub](https://github.com/AparicioJohan/agriutilities) with:

``` r
# install.packages("devtools")
devtools::install_github("AparicioJohan/agriutilities")
```

## Example

This is a basic example which shows you how to use some of the functions
of the package.

The function `check_design_MET` helps us to check the quality of the
data and also to identify the experimental design of the trials. This
works as a quality check or quality control before we fit any model.

``` r
library(agriutilities)
library(agridat)
data(besag.met)
dat <- besag.met
results <- check_design_MET(
  data = dat,
  genotype = "gen",
  trial = "county",
  traits = "yield",
  rep = "rep",
  block = "block",
  col = "col",
  row = "row"
)
```

The results of the previous function are used in `single_model_analysis`
to fit single trial models.

``` r
obj <- single_model_analysis(results, progress = FALSE)
obj$resum_fitted_model
#> $yield
#>   trial heritability        CV    VarGen    VarErr      design
#> 1    C1         0.72  6.057373  85.04747  84.80063 res_row_col
#> 2    C2         0.34 17.906699  23.36095 110.89104 res_row_col
#> 3    C3         0.63 12.189981  82.52294 119.34580 res_row_col
#> 4    C4         0.38  8.222134  32.54677 142.90594 res_row_col
#> 5    C5         0.80  6.805292 104.31185  62.88350 res_row_col
#> 6    C6         0.51 15.267974  71.88271 178.58553 res_row_col
```

The returning object is a set of lists with trial summary, BLUEs, BLUPs,
heritability, variance components, potential extreme observations, and a
lot of more things.
