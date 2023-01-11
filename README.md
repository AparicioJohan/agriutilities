
<!-- README.md is generated from README.Rmd. Please edit that file -->

# agriutilities <img src="man/figures/logo.png" align="right" width="160px"/></a>

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
#> # A tibble: 6 x 9
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
#> # A tibble: 6 x 9
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

The results of the previous function are used in
`single_trial_analysis()` to fit single trial models. This function can
fit, Completely Randomized Designs (**CRD**), Randomized Complete Block
Designs (**RCBD**), Resolvable Incomplete Block Designs (**res-IBD**),
Non-Resolvable Row-Column Designs (**Row-Col**) and Resolvable
Row-Column Designs (**res-Row-Col**).

> **NOTE**: It fits models based on the randomization detected.

``` r
obj <- single_trial_analysis(results, progress = FALSE)
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
#> Online License checked out Wed Jan 11 12:30:38 2023
#> Online License checked out Wed Jan 11 12:30:38 2023
print(met_results)
#> ---------------------------------------------------------------------
#> Trial Effects (BLUEs):
#> ---------------------------------------------------------------------
#>   trait trial predicted.value std.error    status
#> 1 yield    C1       151.40258  1.351961 Estimable
#> 2 yield    C2        65.78758  1.119965 Estimable
#> 3 yield    C3        90.06743  1.478252 Estimable
#> 4 yield    C4       148.09868  1.241288 Estimable
#> 5 yield    C5       122.66271  1.429107 Estimable
#> 6 yield    C6        86.80613  1.492154 Estimable
#> 
#> ---------------------------------------------------------------------
#> Heritability:
#> ---------------------------------------------------------------------
#>   trait        h2
#> 1 yield 0.8244552
#> 
#> ---------------------------------------------------------------------
#> First Overall Predicted Values and Standard Errors (BLUPs):
#> ---------------------------------------------------------------------
#>   trait genotype predicted.value std.error    status
#> 1 yield      G01        110.0200  2.530462 Estimable
#> 2 yield      G02        111.1065  2.541857 Estimable
#> 3 yield      G03        102.7181  2.517908 Estimable
#> 4 yield      G04        115.5364  2.537029 Estimable
#> 5 yield      G05        120.2974  2.544593 Estimable
#> 6 yield      G06        109.4152  2.546403 Estimable
#> 
#> ---------------------------------------------------------------------
#> Variance-Covariance Matrix:
#> ---------------------------------------------------------------------
#> 
#> Correlation Matrix ('us'): yield
#>      C1   C2   C3   C4   C5   C6
#> C1 1.00 0.65 0.58 0.64 0.94 0.40
#> C2 0.65 1.00 0.57 0.79 0.59 0.74
#> C3 0.58 0.57 1.00 0.87 0.74 0.33
#> C4 0.64 0.79 0.87 1.00 0.66 0.41
#> C5 0.94 0.59 0.74 0.66 1.00 0.31
#> C6 0.40 0.74 0.33 0.41 0.31 1.00
#> 
#> Covariance Matrix ('us'): yield
#>       C1    C2    C3    C4     C5    C6
#> C1 74.02 27.52 46.28 32.19  81.34 28.63
#> C2 27.52 24.49 26.31 22.83  29.21 30.15
#> C3 46.28 26.31 86.39 47.71  69.12 25.64
#> C4 32.19 22.83 47.71 34.48  39.07 19.73
#> C5 81.34 29.21 69.12 39.07 101.57 25.94
#> C6 28.63 30.15 25.64 19.73  25.94 68.68
#> 
#> ---------------------------------------------------------------------
#> First Stability Coefficients:
#> ---------------------------------------------------------------------
#>   trait genotype superiority   static    wricke predicted.value
#> 1 yield      G57    22.49186 34.03947 15.764319        92.48602
#> 2 yield      G29    17.19476 34.67881  5.207048        99.16798
#> 3 yield      G59    17.07438 35.32410  4.801998        99.29619
#> 4 yield      G34    16.69884 34.35385  8.138613       100.06035
#> 5 yield      G10    15.76951 33.54824 11.883009       102.00207
#> 6 yield      G31    15.16212 32.21745 10.316367       102.40940
```
