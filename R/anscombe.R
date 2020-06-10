#' Anscombe's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Anscombe61;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model, with or without the
#'    studentising modification of
#'    \insertCite{Bickel78;textual}{skedastic}.
#'
#' Anscombe's Test is among the earliest suggestions for heteroskedasticity
#'    diagnostics in the linear regression model. The test is not based on
#'    formally derived theory but on a test statistic that Anscombe intuited
#'    to be approximately standard normal under the null hypothesis of
#'    homoskedasticity. \insertCite{Bickel78;textual}{skedastic} discusses
#'    the test and suggests a studentising modification (included in this
#'    function) as well as a robustifying modification
#'    (included in \code{\link{bickel}}). The test is two-tailed.
#'
#' The original form of the test is as follows.
#'
#' @param studentise A logical. Should studentising modification of
#'    \insertCite{Bickel78;textual}{skedastic} be implemented? Defaults to
#'    \code{TRUE}; if \code{FALSE}, the original form of the test proposed by
#'    \insertCite{Anscombe61;textual}{skedastic} is used.
#'
#' @inheritParams breusch_pagan
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link{bickel}}, which is a robust extension of this test.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' anscombe(mtcars_lm)
#' anscombe(mtcars_lm, studentise = FALSE)
#'

anscombe <- function (mainlm, studentise = TRUE) {

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
    y <- stats::model.response(stats::model.frame(mainlm))
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

  p <- length(mainlm$coefficients)
  n <- nrow(X)
  yhat <- mainlm$fitted.values
  e <- mainlm$residuals
  M <- fastM(X, n)

  s_sq <- sum(e ^ 2) / (n - p)
  tbar <- sum(diag(M) * yhat) / (n - p)
  if (studentise) {
    method <- "Anscombe-Bickel (studentised)"
    sigma_tilde_sq_xsq <- sum((yhat - tbar) ^ 2) * sum((e ^ 2 -
                          mean(e ^ 2)) ^ 2) / (n - p)
    teststat <- sum(e ^ 2 * (yhat - tbar)) / sqrt(sigma_tilde_sq_xsq)
  } else {
    method <- "Anscombe (non-studentised)"
    sigma_tilde_sq <- 2 * (n - p) / (n - p + 2) * s_sq ^ 2 *
      t(yhat - tbar) %*% (M ^ 2) %*% (yhat - tbar)
    teststat <- sum(e ^ 2 * (yhat - tbar)) / sqrt(sigma_tilde_sq)
  }

  pval <- 2 * stats::pnorm(abs(teststat), lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, parameter = NULL,
               p.value = pval,
               null.value = 0,
               alternative = "two.sided", method = method), class = "htest")
  broom::tidy(rval)
}
