% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/white.R
\name{white}
\alias{white}
\title{White's Test for Heteroskedasticity in a Linear Regression Model}
\usage{
white(mainlm, interactions = FALSE, statonly = FALSE)
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

\item{interactions}{A logical. Should two-way interactions between explanatory
variables be included in the auxiliary regression? Defaults to
\code{FALSE}, since when interaction terms are present the test is not a
pure test of heteroskedasticity but also of model specification.}

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
This function implements the popular method of
   \insertCite{White80;textual}{skedastic} for testing for heteroskedasticity
   in a linear regression model.
}
\details{
White's Test entails fitting an auxiliary regression model in which
   the response variable is the vector of squared residuals from the
   original model and the design matrix includes the original explanatory
   variables, their squares, and (optionally) their two-way interactions.
   The test statistic is the number of observations multiplied by the
   coefficient of determination from the auxiliary regression model:
   \deqn{T = n r_{\mathrm{aux}}^2}

White's Test is thus a special case of the method of
   \insertCite{Breusch79;textual}{skedastic}. Under the null hypothesis of
   homoskedasticity, the distribution of the test statistic is
   asymptotically chi-squared with \code{parameter} degrees of freedom.
   The test is right-tailed.
}
\examples{
mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
white(mtcars_lm)
white(mtcars_lm, interactions = TRUE)

}
\references{
{\insertAllCited{}}
}
\seealso{
This function should not be confused with
  \code{\link[tseries:white.test]{tseries::white.test}}, which does
  \emph{not} implement the method of
  \insertCite{White80;textual}{skedastic} for testing for
  heteroskedasticity in a linear model.
}
