
<!-- README.md is generated from README.Rmd. Please edit that file -->

# agriutilities <a href="https://apariciojohan.github.io/agriutilities/"><img src="man/figures/logo.png" align="right" width="120"/></a>

<!-- badges: start -->
<!-- badges: end -->

The goal of this package is to provide some additional functions for the
analysis of field trials in agriculture. It integrates functionalities
of packages such as `SpATS`, `lme4`, `ASReml-R`, and `statgenSTA`, to
analyse multi-enviromental/trait data in a semi-automatic way.

## Installation

You can install the development version of agriutilities from
[GitHub](https://github.com/AparicioJohan/agriutilities) with:

``` r
# install.packages("devtools")
devtools::install_github("AparicioJohan/agriutilities")
```

## Automatic Data Analysis Pipeline

This is a basic example which shows you how to use some of the functions
of the package.

### Identify the Experimental Design

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

Some of the outputs are:

``` r
print(results)
#> ---------------------------------------------------------------------
#> Summary Traits by Trial:
#> ---------------------------------------------------------------------
#> # A tibble: 6 × 9
#>   county traits  Mean Median    SD    CV     n n_miss miss_perc
#>   <fct>  <chr>  <dbl>  <dbl> <dbl> <dbl> <int>  <int>     <dbl>
#> 1 C1     yield  149.   151.   17.7 0.119   198      6    0.0303
#> 2 C2     yield   56.1   52.1  18.4 0.328   198      6    0.0303
#> 3 C3     yield   87.9   89.2  19.7 0.225   198      6    0.0303
#> 4 C4     yield  145.   143.   17.1 0.118   198      6    0.0303
#> 5 C5     yield  115.   116.   16.4 0.142   198      6    0.0303
#> 6 C6     yield   87.6   87.8  26.6 0.304   198      6    0.0303
#> 
#> ---------------------------------------------------------------------
#> Experimental Design Detected:
#> ---------------------------------------------------------------------
#>   county  exp_design
#> 1     C1 res_row_col
#> 2     C2 res_row_col
#> 3     C3 res_row_col
#> 4     C4 res_row_col
#> 5     C5 res_row_col
#> 6     C6 res_row_col
#> 
#> ---------------------------------------------------------------------
#> Summary Experimental Design:
#> ---------------------------------------------------------------------
#> # A tibble: 6 × 9
#>   county     n n_gen n_rep n_block n_col n_row num_of_reps num_of_gen
#>   <fct>  <int> <int> <int>   <int> <int> <int> <fct>       <fct>     
#> 1 C1       198    64     3       8    11    18 3_9         63_1      
#> 2 C2       198    64     3       8    11    18 3_9         63_1      
#> 3 C3       198    64     3       8    11    18 3_9         63_1      
#> 4 C4       198    64     3       8    11    18 3_9         63_1      
#> 5 C5       198    64     3       8    11    18 3_9         63_1      
#> 6 C6       198    64     3       8    11    18 3_9         63_1      
#> 
#> ---------------------------------------------------------------------
#> Connectivity Matrix:
#> ---------------------------------------------------------------------
#>    C1 C2 C3 C4 C5 C6
#> C1 64 64 64 64 64 64
#> C2 64 64 64 64 64 64
#> C3 64 64 64 64 64 64
#> C4 64 64 64 64 64 64
#> C5 64 64 64 64 64 64
#> C6 64 64 64 64 64 64
```

### Single Trial Analysis (STA)

The results of the previous function are used in `single_model_analysis`
to fit single trial models. (It fits models based on the randomization
detected)

This function can fit, Completely Randomized Designs (**CRD**),
Randomized Complete Block Designs (**RCBD**), Resolvable Incomplete
Block Designs (**res-IBD**), Non-Resolvable Row-Column Designs
(**Row-Col**) and Resolvable Row-Column Designs (**res-Row-Col**)

``` r
obj <- single_model_analysis(results, progress = FALSE)
print(obj)
#> ---------------------------------------------------------------------
#> Summary Fitted Models:
#> ---------------------------------------------------------------------
#>   trait trial heritability        CV    VarGen    VarErr      design
#> 1 yield    C1         0.72  6.057373  85.04747  84.80063 res_row_col
#> 2 yield    C2         0.34 17.906699  23.36095 110.89104 res_row_col
#> 3 yield    C3         0.63 12.189981  82.52294 119.34580 res_row_col
#> 4 yield    C4         0.38  8.222134  32.54677 142.90594 res_row_col
#> 5 yield    C5         0.80  6.805292 104.31185  62.88350 res_row_col
#> 6 yield    C6         0.51 15.267974  71.88271 178.58553 res_row_col
#> 
#> ---------------------------------------------------------------------
#> Outliers Removed:
#> ---------------------------------------------------------------------
#>   trait trial genotype   id outlier
#> 1 yield    C1      G60   50    TRUE
#> 2 yield    C3      G16  534    TRUE
#> 3 yield    C6      G08 1013    TRUE
#> 
#> ---------------------------------------------------------------------
#> First Predicted Values and Standard Errors (BLUEs/BLUPs):
#> ---------------------------------------------------------------------
#>   trait genotype trial    BLUEs  seBLUEs    BLUPs  seBLUPs         wt
#> 1 yield      G01    C1 142.0021 6.678480 144.3256 5.720102 0.02242047
#> 2 yield      G02    C1 160.7412 6.334549 157.4663 5.496182 0.02492118
#> 3 yield      G03    C1 130.8027 6.412041 135.3259 5.558685 0.02432245
#> 4 yield      G04    C1 155.7890 6.459256 154.2345 5.557676 0.02396818
#> 5 yield      G05    C1 163.7357 6.666347 161.2359 5.714016 0.02250216
#> 6 yield      G06    C1 128.2255 6.492415 133.6129 5.608807 0.02372398
```

The returning object is a set of lists with trial summary, BLUEs, BLUPs,
heritability, variance components, potential extreme observations, and a
lot of more things.

### Two-Stage Analysis (MET)

``` r
met_results <- met_analysis(obj)
#> Online License checked out Mon Jan  9 19:29:19 2023
head(met_results$overall_BLUPs)
#>   trait genotype predicted.value std.error    status
#> 1 yield      G01        110.0104  2.531191 Estimable
#> 2 yield      G02        111.0692  2.542598 Estimable
#> 3 yield      G03        102.7300  2.518689 Estimable
#> 4 yield      G04        115.5712  2.537870 Estimable
#> 5 yield      G05        120.2970  2.544906 Estimable
#> 6 yield      G06        109.4588  2.546752 Estimable
```
