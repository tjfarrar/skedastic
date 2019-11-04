# skedastic

The purpose of the `skedastic` package is to make accessible diagnostic methods 
for detecting heteroskedasticity in linear regression models.

## Installation

```R
# Install from CRAN
install.packages("skedastic", dependencies = c("Depends", "Imports"))

# Or the development version from GitHub:
install.packages("devtools")
devtools::install_github("tjfarrar/skedastic")
```

## Usage

The purpose of the `skedastic` package is to make accessible diagnostic methods 
for detecting heteroskedasticity in linear regression models. Heteroskedasticity 
(sometimes spelt 'heteroscedasticity') is a violation of one of the assumptions 
of the classical linear regression model (the Gauss-Markov Assumptions). This assumption, known as homoskedasticity, holds that the variance of the random 
error term remains constant across all observations.

Most of the functions in the package implement particular hypothesis testing 
methods for detecting heteroskedasticity that have been published in a journal 
article or book. Other functions implement graphical methods for detecting 
heteroskedasticity or compute transformations of the Ordinary Least Squares 
(OLS) residuals that are useful in heteroskedasticity detection.

Nearly all of the exported functions in the package require as an argument an  
`lm` object (or, for advanced users desiring faster computation, a list 
containing the response vector and the design matrix in the exact form 
required by `lm.fit`). Thus, to use this package a user must first be familiar 
with how to fit linear regression models using the `lm` function from package 
`stats`. Note that the package currently supports only linear regression models 
fit using OLS.

For heteroskedasticity tests that are implemented in other R packages on CRAN, 
or in other statistical software such as SAS or SHAZAM, the functions in the 
`skedastic` package have been checked against them to ensure that they produce 
the same values of the test statistic and $p$-value. This is true for the 
following functions: `breusch_pagan`, `cook_weisberg`, `glejser`, 
`goldfeld_quandt` (parametric test only), `harvey`, and `white_lm`.

Here is an example of implementing the Breusch-Pagan Test for heteroskedasticity 
on a linear regression model fit to the [`cars`](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/cars.html) 
dataset, with distance (`cars$dist`) as the response (dependent) variable and 
speed (`cars$speed`) as the explanatory (independent) variable.

```R
library(skedastic)
mylm <- lm(dist ~ speed, data = cars)
bp <- breusch_pagan(mylm)
bp
attributes(bp)
bp$statistic
bp$p.value
```

To compute BLUS residuals for the same model:

```R
myblusres <- blus(mylm, omit = "last")
myblusres
# Plot OLS residuals against BLUS residuals
plot(mylm$residuals, myblusres)
```

To create customised residual plots for the same model:

```R
hetplot(mylm, horzvar = c("explanatory", "log_explanatory"), vertvar = c("res", "res_stud"), vertfun = "2", filetype = NULL)
```

## Learn More

No vignettes have been created yet for this package. Watch this space.
