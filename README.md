
<!-- README.md is generated from README.Rmd. Please edit that file -->

# playbase

<!-- badges: start -->
<!-- badges: end -->

The \`playbase´ package contains the core back-end functionality for the
OmicsPlayground. This package allows you to run, develop, and test any
essential functions used in the OmicsPlayground directly in the R
console without needing to worry about the R Shiny front-end code.

## Installation

You can install the development version of playbase from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bigomics/playbase")
```

## Data upload

The first step in any OmicsPlayground analysis is to upload data which
can be used to create a `pgx` object. The `pgx` object is basically the
core data structure in the OmicsPlayground upon which most analysis and
plotting functions operate.

``` r
library(playbase)

## example files can be accessed via playbase_example()
counts <- playbase::read_counts(playbase::example_file('counts.csv'))
#> 
#> 
samples <- playbase::read_samples(playbase::example_file('samples.csv'))
contrasts <- playbase::read_contrasts(playbase::example_file('contrasts.csv'))
```
