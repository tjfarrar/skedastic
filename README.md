# skedastic

[![Build Status](https://travis-ci.org/tjfarrar/skedastic.svg?branch=master)](https://travis-ci.org/tjfarrar/skedastic)

The purpose of the `skedastic` package is to make diagnostic methods 
for detecting heteroskedasticity in linear regression models accessible to R 
users.

## Installation

```R
# Install from CRAN
install.packages("skedastic", dependencies = c("Depends", "Imports"))

# Or the development version from GitHub:
install.packages("devtools")
devtools::install_github("tjfarrar/skedastic")
```

## Usage

The purpose of the `skedastic` package is to make diagnostic methods 
for detecting heteroskedasticity in linear regression models accessible to R 
users. Heteroskedasticity (sometimes spelt 'heteroscedasticity') is a violation 
of one of the assumptions of the classical linear regression model (the 
Gauss-Markov Assumptions). This assumption, known as homoskedasticity, holds 
that the variance of the random error term remains constant across all 
observations.

23 distinct functions in the package implement hypothesis testing 
methods for detecting heteroskedasticity that have been previously published in academic literature. Other functions implement graphical methods for detecting 
heteroskedasticity or perform supporting tasks for the tests such as computing transformations of the Ordinary Least Squares (OLS) residuals that are useful in heteroskedasticity detection, or computing probabilities from the null distribution of a nonparametric test statistic. Certain functions have applications beyond the problem of heteroskedasticity in linear regression. These include `pRQF`, which computes cumulative probabilities from the distribution of a ratio of quadratic forms in normal random vectors, `twosidedpval`, which implements three different approaches for calculating two-sided $p$-values from asymmetric null distributions, and `dDtrend` and `pdDtrend`, which compute probabilities from Lehmann's nonparametric trend statistic. 

Most of the exported functions in the package take a linear model as their primary argument (which can be passed as an `lm` object). Thus, to use this package a user must first be familiar 
with how to fit linear regression models using the `lm` function from package 
`stats`. Note that the package currently supports only linear regression models 
fit using OLS.

For heteroskedasticity tests that are implemented in other R packages on CRAN, 
or in other statistical software such as SAS or SHAZAM, the functions in the 
`skedastic` package have been checked against them to ensure that they produce 
the same values of the test statistic and $p$-value. This is true of `breusch_pagan`, `cook_weisberg`, `glejser`, 
`goldfeld_quandt` (parametric test only), `harvey`, and `white_lm`.

Here is an example of implementing the Breusch-Pagan Test for heteroskedasticity 
on a linear regression model fit to the [`cars`](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/cars.html) 
dataset, with distance (`cars$dist`) as the response (dependent) variable and 
speed (`cars$speed`) as the explanatory (independent) variable.

```R
library(skedastic)
mylm <- lm(dist ~ speed, data = cars)
breusch_pagan(mylm)
```

To compute BLUS residuals for the same model:

```R
myblusres <- blus(mylm, omit = "last")
myblusres
```

To create customised residual plots for the same model:

```R
hetplot(mylm, horzvar = c("explanatory", "log_explanatory"), vertvar = c("res", "res_stud"), vertfun = "2", filetype = NA)
```

## Learn More

No vignettes have been created yet for this package. Watch this space.
