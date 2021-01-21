
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: covidprobability

<!-- badges: start -->

<!-- badges: end -->

This package provides the functions, data and documentation that support
a calculator to determine the probability of an undetected COVID-19
infection in a setting/unit after a potential exposure, testing, and
when there are no symptomatic cases. For a detailed explanation of the
rationale and implementation, please see the vignette.

## Installation

Warning this package is still in development.

You can install the latest version of covidprobability from
[Github](https://github.com/eebrown/covidprobability) with:

``` r
devtools::install.github("eebrown/covidprobability")
```

## Example

``` r
library(covidprobability)

test_n <- unit_probability(test_day = 9, pre0 =  0.13, sens = sens, spec = 1, 
                           incubation = incubation, asympt = 0.278, n = 10)
```

<img src="man/figures/README-unit_example-1.png" width="100%" />