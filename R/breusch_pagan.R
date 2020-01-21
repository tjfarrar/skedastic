#' Breusch-Pagan Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the popular method of
#'    \insertCite{Breusch79;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model, with or without the
#'    studentising modification of \insertCite{Koenker81;textual}{skedastic}.
#'
#' The Breusch-Pagan Test entails fitting an auxiliary regression model in
#'    which the response variable is the vector of squared residuals from the
#'    original model and the design matrix \eqn{Z} consists of one or more
#'    exogenous variables that are suspected of being related to the error
#'    variance. In the absence of prior information on a possible choice of
#'    \eqn{Z}, one would typically use the explanatory variables from the
#'    original model. Under the null hypothesis of homoskedasticity, the
#'    distribution of the test statistic is asymptotically chi-squared with
#'    \code{parameter} degrees of freedom. The test is right-tailed.
#' @param mainlm Either an object of \code{\link[base]{class}} "lm", or a list
#'    of two components: a response vector and a design matrix (in that order).
#'    If the latter, the design matrix must begin with a column of 1s if an
#'    intercept is to be included in the linear model. Passing an object of
#'    class "lm" is recommended in applications; passing a list containing
#'    the data is recommended where optimising computational performance is
#'    important.
#' @param auxdesign A \code{\link[base]{data.frame}} or
#'    \code{\link[base]{matrix}} representing an auxiliary design matrix of
#'    containing exogenous variables that (under alternative hypothesis) are
#'    related to error variance, or a character "fitted.values" indicating
#'    that the fitted \eqn{\hat{y}_i} values from OLS should be used.
#'    If set to \code{NULL} (the default), the
#'    design matrix of the original regression model is used. An intercept
#'    is included in the auxiliary regression even if the first column of
#'    \code{auxdesign} is not a vector of ones.
#' @param koenker A logical. Should studentising modification of
#'    \insertCite{Koenker81;textual}{skedastic} be implemented? Defaults to
#'    \code{TRUE}; if \code{FALSE}, the original form of the test proposed by
#'    \insertCite{Breusch79;textual}{skedastic} is used.
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[lmtest:bptest]{lmtest::bptest}}, which performs exactly
#'    the same test as this function; \code{\link[car:ncvTest]{car::ncvTest}},
#'    which is not the same test but is implemented in
#'    \code{\link{cook_weisberg}}; \code{\link{white_lm}}, which is a special
#'    case of the Breusch-Pagan Test.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' breusch_pagan(mtcars_lm)
#' breusch_pagan(mtcars_lm, koenker = FALSE)
#'

breusch_pagan <- function (mainlm, auxdesign = NULL, koenker = TRUE) {

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x), is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, drop = FALSE]
    }
    mainlm <- stats::lm.fit(X, y)
  }

  if (is.null(auxdesign)) {
    Z <- X
  } else if (is.character(auxdesign)) {
    if (auxdesign == "fitted.values") {
      Z <- t(t(mainlm$fitted.values))
    } else stop("Invalid character value for `auxdesign`")
  } else {
    Z <- auxdesign
    if (nrow(auxdesign) != nrow(X)) stop("No. of observations in `auxdesign`
                                         must match\nno. of observations in
                                         original model.")
  }

  hasintercept <- columnof1s(Z)
  if (!hasintercept[[1]]) {
    Z <- cbind(1, Z)
    message("Column of 1's added to `auxdesign`")
  }

  p <- ncol(Z) - 1
  n <- nrow(Z)
  sigma_hatsq <- sum(mainlm$residuals ^ 2) / n
  w_hat <- mainlm$residuals ^ 2 - sigma_hatsq

  if (koenker) {
    method <- "Koenker (studentised)"
    teststat <- n * sum(stats::lm.fit(Z, w_hat)$fitted.values ^ 2) / sum(w_hat ^ 2)
  } else {
    method <- "Breusch-Pagan (non-studentised)"
    teststat <- sum(stats::lm.fit(Z, w_hat)$fitted.values ^ 2) / (2 * sigma_hatsq ^ 2)
  }
  pval <- 1 - stats::pchisq(teststat, df = p)

  rval <- structure(list(statistic = teststat, parameter = p, p.value = pval,
               null.value = "Homoskedasticity",
               alternative = "Heteroskedasticity", method = method), class = "htest")
  broom::tidy(rval)
}
