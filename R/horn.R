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
#'    improves the normal approximation. When
#'    exact distribution is used, two-sided \eqn{p}-value is computed by
#'    doubling the one-sided \eqn{p}-value, since the distribution of \eqn{D}
#'    is symmetric. The function does not support the exact distribution of
#'    \eqn{D} in the presence of ties, so in this case the normal approximation
#'    is used regardless of \eqn{n}.
#'
#' @param restype A character specifying which residuals to use: OLS residuals
#'    (the default) or the \link[=blus]{BLUS residuals} of
#'    \insertCite{Theil65;textual}{skedastic}.
#' @param alternative A character specifying the form of alternative
#'    hypothesis; one of "two.sided" (default), "greater", or "less".
#'    "two.sided" corresponds to any trend in the absolute residuals when
#'    ordered by \code{deflator}. "greater" corresponds to a negative trend in
#'    the absolute residuals when ordered by \code{deflator}. "less"
#'    corresponds to a positive trend in the absolute residuals when ordered by
#'    \code{deflator}. (Notice that \eqn{D} tends to be small when there is a
#'    positive trend.)
#' @param exact A logical. Should exact \eqn{p}-values be computed? If
#'    \code{FALSE} (the default), the normalised statistic \eqn{Z}
#' @param correct A logical. Should a continuity correction be used when
#'    computing the \eqn{p}-value? This parameter is ignored if
#'    \code{exact == TRUE} or if ties are present in the ranks.
#' @param ... Optional further arguments to pass to \code{\link{blus}}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams goldfeld_quandt
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
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

horn <- function (mainlm, deflator = NULL, restype = c("ols", "blus"),
                  alternative = c("two.sided", "greater", "less"),
                  exact = (n <= 10), correct = !exact, ...) {

  restype <- match.arg(restype, c("ols", "blus"))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
    p <- ncol(X)
    hasintercept <- columnof1s(X)
    y <- stats::model.response(stats::model.frame(mainlm))
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x),
                                        is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, ]
    }
    p <- ncol(X)
    hasintercept <- columnof1s(X)
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)

  if (is.numeric(deflator) && deflator == as.integer(deflator)) {
    if (hasintercept[[1]] && deflator == 1) {
      stop("deflator cannot be the model intercept")
    } else if (deflator > p) {
      stop("`deflator` is not the index of a column of design matrix")
    }
  } else if (is.character(deflator)) {
    if (deflator == "(Intercept)") {
      stop("deflator cannot be the model intercept")
    } else if (!deflator %in% colnames(X)) {
      stop("`deflator` is not the name of a column of design matrix")
    }
  } else if (!is.null(deflator)) stop("`deflator` must be integer or character")

  if (is.null(deflator)) {
    if (restype == "ols") {
      absres <- abs(mainlm$residuals)
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...))
      absres <- absres[!is.na(absres)]
      n <- length(absres)
    }
  } else {
    if (restype == "ols") {
      absres <- abs(mainlm$residuals)[order(X[, deflator])]
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...)[order(X[, deflator])])
      absres <- absres[!is.na(absres)]
      n <- length(absres)
    }
  }

  R <- rank(absres, ties.method = "average")
  teststat <- sum((R - 1:n) ^ 2)
  ties <- ifelse(sum(R) != n * (n + 1) / 2, TRUE, FALSE)
  if (ties) {
    if (exact) warning("Ties are present and exact distribution is not
                       available for D in presence of ties. Normal
                       approximation will be used.")

    d <- as.double(names(table(R)))
    if (max(d) / n == 1) warning("Normal approximation may not be accurate
                                 since maximum rank / n = 1. See
                                 Lehmann (1975), p. 294.")
    ED <- (n ^ 3 - n) / 6 - sum(d ^ 3 - d) / 12
    VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36 * (1 - sum(d ^ 3 - d) /
                                                    (n ^ 3 - n))
    teststat <- (teststat - ED) / sqrt(VD)
    pval <- switch(alternative, "greater" = stats::pnorm(teststat,
                  lower.tail = FALSE), "less" = stats::pnorm(teststat,
                  lower.tail = TRUE), "two.sided" =
                     2 * stats::pnorm(abs(teststat),
                                      lower.tail = FALSE))
  } else {
    pval <- switch(alternative, "greater" = pDtrend(n, teststat,
                                      lower.tail = FALSE, exact, correct),
                                "less" = pDtrend(n, teststat,
                                      lower.tail = TRUE, exact, correct),
                                "two.sided" = 2 * pDtrend(n, abs(teststat),
                                      lower.tail = ifelse(exact,
                                        teststat < (n * (n - 1) * (n + 1) / 6),
                                        FALSE), exact, correct))
  }

  rval <- structure(list(statistic = teststat, p.value = pval, parameter = n,
               null.value = "Homoskedasticity", alternative = alternative),
               class = "htest")
  broom::tidy(rval)
}
