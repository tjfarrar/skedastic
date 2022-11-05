# skedastic

[![Build Status](https://travis-ci.org/tjfarrar/skedastic.svg?branch=master)](https://travis-ci.org/tjfarrar/skedastic)

The purpose of the `skedastic` package is to make a suite of old and new 
methods for detecting and correcting for heteroskedasticity in linear 
regression models accessible to R users.

## Installation

```R
# Install from CRAN
install.packages("skedastic", dependencies = c("Depends", "Imports"))

# Or the development version from GitHub:
install.packages("devtools")
devtools::install_github("tjfarrar/skedastic")
```

## Usage

Heteroskedasticity (sometimes spelt 'heteroscedasticity') is a violation 
of one of the assumptions of the classical linear regression model (the 
Gauss-Markov Assumptions). This assumption, known as homoskedasticity, holds 
that the variance of the random error term remains constant across all 
observations. Under heteroskedasticity, the Ordinary Least Squares estimator 
is no longer the Best Linear Unbiased Estimator (BLUE) of the parameter 
vector, while the classical t-tests for testing significance of the parameters 
are invalid. Thus, heteroskedasticity-robust methods are required.

The most novel functionality of this package is provided by the `alvm.fit` and 
`anlvm.fit` functions, which fit an Auxiliary Linear Variance Model or an 
Auxiliary Nonlinear Variance Model, respectively. These are new models for 
estimating error variances in heteroskedastic linear regression models, 
developed as part of the author's doctoral research.

The `hccme` function computes heteroskedasticity-consistent covariance matrix 
estimates for the $\hat{\beta}$ Ordinary Least Squares estimator using ten 
different methods found in the literature.

25 distinct functions in the package implement hypothesis testing 
methods for detecting heteroskedasticity that have been previously published in academic literature. Other functions implement graphical methods for detecting 
heteroskedasticity or perform supporting tasks for the tests such as computing transformations of the Ordinary Least Squares (OLS) residuals that are useful in heteroskedasticity detection, or computing probabilities from the null distribution of a nonparametric test statistic. Certain functions have applications beyond the problem of heteroskedasticity in linear regression. These include `pRQF`, which computes cumulative probabilities from the distribution of a ratio of quadratic forms in normal random vectors, `twosidedpval`, which implements three different approaches for calculating two-sided $p$-values from asymmetric null distributions, and `dDtrend` and `pdDtrend`, which compute probabilities from Lehmann's nonparametric trend statistic. 

Most of the exported functions in the package take a linear model as their primary argument (which can be passed as an `lm` object). Thus, to use this package a user must first be familiar with how to fit linear regression models using the `lm` function from package `stats`.

Here is an example of implementing the Breusch-Pagan Test for heteroskedasticity on a linear regression model fit to the [`cars`](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/cars.html) 
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

To fit an auxiliary linear variance model to the same linear regression model,  assuming that the error variances are a linear function of the `speed` predictor, and extract the resulting variance estimates:

```R
myalvm <- alvm.fit(mylm, model = "linear")
myalvm$var.est
```

To fit an auxiliary linear variance model to the same linear regression model,  assuming that the error variances are a quadratic function of the `speed` predictor, and extract the resulting variance estimates:

```R
myanlvm <- anlvm.fit(mylm, g = function(x) x ^ 2)
mynalvm$var.est
```

## Learn More

No vignettes have been created yet for this package. Watch this space.
