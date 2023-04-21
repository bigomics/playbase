
<!-- README.md is generated from README.Rmd. Please edit that file -->

# playbase

<!-- badges: start -->

[![R package
check](https://github.com/bigomics/playbase/actions/workflows/r.yml/badge.svg)](https://github.com/bigomics/playbase/actions/workflows/r.yml)
[![Codecov test
coverage](https://codecov.io/gh/bigomics/playbase/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bigomics/playbase?branch=master)

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
#> 
#> 
#> 

## example files can be accessed via the example_file() function
counts <- playbase::read_counts(playbase::example_file("counts.csv"))
samples <- playbase::read_samples(playbase::example_file("samples.csv"))
contrasts <- playbase::read_contrasts(playbase::example_file("contrasts.csv"))

## create a pgx object (new style)
#pgx <- playbase::create_pgx(counts, samples, contrasts)

## old style still works
#pgx <- playbase::pgx.createPGX(counts, samples, contrasts)
```

## Adding namespace

``` r
library(addns)
#addns::add_ns(path = "R/", packages = c("playbase"))
#adds::rm_imports(path = "R/")
```

## The pgx object

The core object in playbase is the `pgx` object. This object holds the
raw data and any analysis results returned from playbase modules /
boards. The pgx object is simply an R list. It contains minimally the
following list items:

- counts
- samples
- contrasts

A pgx object is created from these three list items via the following
function:

    my_pgx <- pgx.createPGX(counts, samples, contrasts)

Once a pgx object is created from these three items, the various
playbase modules can operate on the pgx object to generate the analysis
results relevant to that specific module.

## Playbase modules

As mentioned above, the core object in playbase is the `pgx` object.
This holds all of the analysis and results derived from the raw data, as
well as the raw data itself. There are various modules in playbase that
take a pgx object as input, perform some analysis on the raw data in the
pgx object, and then append these results to the pgx object. These
modules are more-or-less independent of one another and can therefore be
parallelized or run in any arbitrary order.

The core playbase modules operate on either genes or genesets.

The gene modules are as follows:

- ttest
- ttest.welch
- ttest.rank
- voom.limma
- trend.limma
- notrend.limma
- edger.qlf
- edger.lrt
- deseq2.wald
- deseq2.lrt

The geneset methods are as follows:

- fisher
- gsva
- ssgsea
- spearman
- camera
- fry
- fgsea

And extra modules are as follows:

- meta.go
- deconv
- infer
- drugs
- connectivity
- wordcloud
