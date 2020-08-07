#' Horn's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the nonparametric test of
#'    \insertCite{Horn81;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The test entails specifying a 'deflator', an explanatory variable
#'    suspected of being related to the error variance. Residuals are ordered
#'    by the deflator and the nonparametric trend statistic
#'    \eqn{D=\sum (R_i - i)^2} proposed by
#'    \insertCite{Lehmann75;textual}{skedastic} is
#'    then computed on the absolute residuals and used to test for an
#'    increasing or decreasing trend, either of which would correspond to
#'    heteroskedasticity. Exact probabilities for the null distribution of
#'    \eqn{D} can be obtained from functions \code{\link{dDtrend}} and
#'    \code{\link{pDtrend}}, but since computation time increases rapidly with
#'    \eqn{n}, use of a normal approximation is recommended for \eqn{n>10}.
#'    \insertCite{Lehmann75;textual}{skedastic} proves that \eqn{D} is
#'    asymptotically normally distributed and the approximation of the
#'    statistic \eqn{Z=(D-E(D))/\sqrt{V(D)}} to the standard normal
#'    distribution is already quite good for \eqn{n=11}.
#'
#'    The expectation and variance of \eqn{D} (when ties are absent) are
#'    respectively \eqn{E(D)=\frac{n^3-n}{6}} and
#'    \eqn{V(D)=\frac{n^2(n+1)^2(n-1)}{36}}; see
#'    \insertCite{Lehmann75;textual}{skedastic} for \eqn{E(D)} and \eqn{V(D)}
#'    when ties are present. When ties are absent, a continuity correction
#'    is used to improve the normal approximation. When
#'    exact distribution is used, two-sided \eqn{p}-value is computed by
#'    doubling the one-sided \eqn{p}-value, since the distribution of \eqn{D}
#'    is symmetric. The function does not support the exact distribution of
#'    \eqn{D} in the presence of ties, so in this case the normal approximation
#'    is used regardless of \eqn{n}.
#'
#' @param alternative A character specifying the form of alternative
#'    hypothesis; one of \code{"two.sided"} (default), \code{"greater"},
#'    or \code{"less"}. \code{"two.sided"} corresponds to any trend in the
#'    absolute residuals when ordered by \code{deflator}. \code{"greater"}
#'    corresponds to a negative trend in the absolute residuals when ordered by
#'    \code{deflator}. \code{"less"} corresponds to a positive trend in the
#'    absolute residuals when ordered by \code{deflator}. (Notice that \eqn{D}
#'    tends to be small when there is a positive trend.)
#' @param exact A logical. Should exact \eqn{p}-values be computed? If
#'    \code{FALSE}, a normal approximation is used. Defaults to \code{TRUE}
#'    only if the number of absolute residuals being ranked is \eqn{\le 10}.
#' @param restype A character specifying which residuals to use: \code{"ols"}
#'    for OLS residuals (the default) or the \code{"blus"} for
#'    \link[=blus]{BLUS} residuals. The advantage of using BLUS residuals is
#'    that, under the null hypothesis, the assumption that the random series
#'    is independent and identically distributed is met (whereas with OLS
#'    residuals it is not). The disadvantage of using BLUS residuals is that
#'    only \eqn{n-p} residuals are used rather than the full \eqn{n}.
#' @param ... Optional further arguments to pass to \code{\link{blus}}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams goldfeld_quandt
#'
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' horn(mtcars_lm, deflator = "qsec")
#' horn(mtcars_lm, deflator = "qsec", restype = "blus")
#'

horn <- function(mainlm, deflator = NA, restype = c("ols", "blus"),
                  alternative = c("two.sided", "greater", "less"),
                  exact = (nres <= 10), statonly = FALSE, ...) {

  restype <- match.arg(restype, c("ols", "blus"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  processmainlm(m = mainlm, needy = FALSE)

  hasintercept <- columnof1s(X)
  if (class(mainlm) == "list") {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)

  checkdeflator(deflator, X, p, hasintercept[[1]])

  if (is.na(deflator) || is.null(deflator)) {
    if (restype == "ols") {
      absres <- abs(e)
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...))
      absres <- absres[!is.na(absres)]
    }
  } else {
    if (!is.na(suppressWarnings(as.integer(deflator)))) {
      deflator <- as.integer(deflator)
    }
    if (restype == "ols") {
      absres <- abs(e)[order(X[, deflator])]
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...)[order(X[, deflator])])
      absres <- absres[!is.na(absres)]
    }
  }
  nres <- length(absres)

  R <- data.table::frank(absres, ties.method = "average")
  teststat <- sum((R - 1:nres) ^ 2)
  if (statonly) return(teststat)
  d <- table(R)
  is.ties <- (max(d) > 1)

  if (is.ties) {
    if (exact) {
      warning("Ties are present and exact distribution is not available for D in presence of ties. Normal approximation will be used.")
      exact <- FALSE
    }
    if (max(d) / nres == 1) warning("Normal approximation may not be accurate since maximum rank / m = 1. See Lehmann (1975), p. 294.")
  } else {
    d <- NULL
  }

  twosided <- function(tstat, m) {
    if (teststat > (m * (m - 1) * (m + 1) / 6)) {
      return(2 * pDtrend(k = tstat, n = m, lower.tail = FALSE,
                         exact, tiefreq = d))
    } else if (tstat < (m * (m - 1) * (m + 1) / 6)) {
      return(2 * pDtrend(k = tstat, n = m, lower.tail = TRUE,
                         exact, tiefreq = d))
    } else if (tstat == (m * (m - 1) * (m + 1) / 6)) {
      return(1)
    }
  }
  pval <- switch(alternative, "greater" = pDtrend(k = teststat, n = nres,
                 lower.tail = FALSE, exact, tiefreq = d),
                 "less" = pDtrend(k = teststat, n = nres,
                 lower.tail = TRUE, exact, tiefreq = d),
                 "two.sided" = twosided(teststat, m = nres))

  rval <- structure(list(statistic = teststat, p.value = pval,
               alternative = alternative),
               class = "htest")
  broom::tidy(rval)
}
