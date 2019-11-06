#' White's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the popular method of
#'    \insertCite{White80;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' White's Test entails fitting an auxiliary regression model in which the
#' response variable is the vector of squared residuals from the original
#' model and the design matrix includes the original explanatory variables,
#' their squares, and (optionally) their two-way interactions. The test
#' statistic is the number of observations multiplied by the coefficient of
#' determination from the auxiliary regression model:
#' \deqn{T = n r_{\mathrm{aux}}^2}
#' White's Test is thus a special case of the method of
#' \insertCite{Breusch79;textual}{skedastic}. Under the null hypothesis of
#'    homoskedasticity, the distribution of the test statistic is
#'    asymptotically chi-squared with \code{parameter} degrees of freedom.
#'    The test is right-tailed.
#' @param interactions A logical. Should two-way interactions between explanatory
#'    variables be included in the auxiliary regression? Defaults to
#'    \code{FALSE}, since when interaction terms are present the test is not a
#'    pure test of heteroskedasticity but also of model specification.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso This function should not be confused with
#'   \code{\link[het.test:whites.htest]{het.test::whites.htest}} and
#'   \code{\link[tseries:white.test]{tseries::white.test}}, which do \emph{not}
#'   implement the method of \insertCite{White80;textual}{skedastic} for
#'   testing for heteroskedasticity in a linear model.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' white_lm(mtcars_lm)
#' white_lm(mtcars_lm, interactions = TRUE)
#'

white_lm <- function (mainlm, interactions = FALSE) {

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x), is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, ]
    }
    mainlm <- stats::lm.fit(X, y)
  }

  hasintercept <- columnof1s(X)
  if (!hasintercept[[1]]) {
    message("Intercept included in auxiliary design matrix")
  } else {
    X <- X[, -hasintercept[[2]]]
  }

  n <- nrow(X)
  k <- ncol(X)

  Z <- cbind(1, X, X ^ 2)
  if (interactions) {
    Z <- cbind(Z, generate_interactions(X))
  }

  esq <- mainlm$residuals ^ 2
  auxlm <- stats::lm.fit(Z, esq)
  iota <- rep(1, n)
  N <- diag(n) - 1 / n * (iota %*% t(iota))
  e_aux <- auxlm$residuals
  teststat <- n * (1 - (t(e_aux) %*% e_aux) / (t(esq) %*% N %*% esq))

  df <- ncol(Z) - 1
  pval <- 1 - stats::pchisq(teststat, df = df)

  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "Heteroskedasticity",
                         method = "White's Test"), class = "htest")
  broom::tidy(rval)
}
