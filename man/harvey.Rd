% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harvey.R
\name{harvey}
\alias{harvey}
\title{Harvey Test for Heteroskedasticity in a Linear Regression Model}
\usage{
harvey(mainlm, auxdesign = NA, statonly = FALSE)
}
\arguments{
\item{mainlm}{Either an object of \code{\link[base]{class}} \code{"lm"}
(e.g., generated by \code{\link[stats]{lm}}), or
a list of two objects: a response vector and a design matrix. The objects
are assumed to be in that order, unless they are given the names
\code{"X"}  and \code{"y"} to distinguish them. The design matrix passed
in a list must begin with a column of ones if an intercept is to be
included in the linear model. The design matrix passed in a list should
not contain factors, as all columns are treated 'as is'. For tests that
use ordinary least squares residuals, one can also pass a vector of
residuals in the list, which should either be the third object or be
named \code{"e"}.}

\item{auxdesign}{A \code{\link[base]{data.frame}} or
\code{\link[base]{matrix}} representing an auxiliary design matrix of
containing exogenous variables that (under alternative hypothesis) are
related to error variance, or a character "fitted.values" indicating
that the fitted \eqn{\hat{y}_i} values from OLS should be used.
If set to \code{NA} (the default), the
design matrix of the original regression model is used. An intercept
is included in the auxiliary regression even if the first column of
\code{auxdesign} is not a vector of ones.}

\item{statonly}{A logical. If \code{TRUE}, only the test statistic value
is returned, instead of an object of \code{\link[base]{class}}
\code{"htest"}. Defaults to \code{FALSE}.}
}
\value{
An object of \code{\link[base]{class}} \code{"htest"}. If object is
   not assigned, its attributes are displayed in the console as a
   \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
}
\description{
This function implements the method of
   \insertCite{Harvey76;textual}{skedastic} for testing for "multiplicative"
   heteroskedasticity in a linear regression model.
   \insertCite{Mittelhammer00;textual}{skedastic} gives the
   formulation of the test used here.
}
\details{
Harvey's Test entails fitting an auxiliary regression model in
   which the response variable is the log of the  vector of squared
   residuals from the original model and the design matrix \eqn{Z}
   consists of one or more exogenous variables that are suspected of being
   related to the error variance. In the absence of prior information on
   a possible choice of \eqn{Z}, one would typically use the explanatory
   variables from the original model. Under the null hypothesis of
   homoskedasticity, the distribution of the test statistic is
   asymptotically chi-squared with \code{parameter} degrees of freedom.
   The test is right-tailed.
}
\examples{
mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
harvey(mtcars_lm)
harvey(mtcars_lm, auxdesign = "fitted.values")

}
\references{
{\insertAllCited{}}
}
\seealso{
the description of the test in
  \href{http://www.econometrics.com/intro/testhet.htm}{SHAZAM} software
  (which produces identical results).
}
