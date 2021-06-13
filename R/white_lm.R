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
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso This function should not be confused with
#'   \code{\link[tseries:white.test]{tseries::white.test}}, which does
#'   \emph{not} implement the method of
#'   \insertCite{White80;textual}{skedastic} for testing for
#'   heteroskedasticity in a linear model.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' white_lm(mtcars_lm)
#' white_lm(mtcars_lm, interactions = TRUE)
#'

white_lm <- function(mainlm, interactions = FALSE, statonly = FALSE) {

  processmainlm(m = mainlm, needy = FALSE)

  hasintercept <- columnof1s(X)
  if (!hasintercept[[1]]) {
    message("Intercept included in auxiliary design matrix")
  } else {
    X <- X[, -hasintercept[[2]], drop = FALSE]
  }

  n <- nrow(X)

  Z <- cbind(1, X, X ^ 2)
  if (interactions) {
    Z <- cbind(Z, generate_interactions(X))
  }

  esq <- e ^ 2
  auxlm <- stats::lm.fit(Z, esq)
  iota <- rep(1, n)
  N <- diag(n) - 1 / n * (tcrossprod(iota))
  e_aux <- auxlm$residuals
  teststat <- as.double(n * (1 - crossprod(e_aux) / (t(esq) %*% N %*% esq)))
  if (statonly) return(teststat)

  df <- ncol(Z) - 1
  pval <- stats::pchisq(teststat, df = df, lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater",
                         method = "White's Test"), class = "htest")
  broom::tidy(rval)
}
