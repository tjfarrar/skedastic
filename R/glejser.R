#' Glejser Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Glejser69;textual}{skedastic} for testing for "multiplicative"
#'    heteroskedasticity in a linear regression model.
#'    \insertCite{Mittelhammer00;textual}{skedastic} gives the
#'    formulation of the test used here.
#'
#' @details Glejser's Test entails fitting an auxiliary regression model in
#'    which the response variable is the absolute residual from the original
#'    model and the design matrix \eqn{Z} consists of one or more exogenous
#'    variables that are suspected of being related to the error variance.
#'    In the absence of prior information on a possible choice of \eqn{Z},
#'    one would typically use the explanatory variables from the original model.
#'    Under the null hypothesis of homoskedasticity, the distribution of the
#'    test statistic is asymptotically chi-squared with \code{parameter} degrees
#'    of freedom. The test is right-tailed.
#'
#' @param sigmaest A character indicating which model residuals to use in the
#'    \eqn{\hat{\omega}} estimator in the denominator of the test statistic.
#'    If \code{"main"} (the default), the OLS residuals from the original model
#'    are used; this produces results identical to the Glejser Test in SHAZAM
#'    software. If \code{"auxiliary"}, the OLS residuals from the auxiliary
#'    model are used, as in \insertCite{Mittelhammer00;textual}{skedastic}.
#'    Partial matching is used.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso the description of the test in
#'   \href{http://www.econometrics.com/intro/testhet.htm}{SHAZAM} software
#'   (which produces identical results).
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' glejser(mtcars_lm)
#'

glejser <- function(mainlm, auxdesign = NA,
                     sigmaest = c("main", "auxiliary"), statonly = FALSE) {

  sigmaest <- match.arg(sigmaest, c("main", "auxiliary"))

  auxfitvals <- ifelse(all(is.na(auxdesign)) | is.null(auxdesign), FALSE,
                                    auxdesign == "fitted.values")
  processmainlm(m = mainlm, needy = auxfitvals, needyhat = auxfitvals,
                needp = FALSE)

  if (all(is.na(auxdesign)) || is.null(auxdesign)) {
    Z <- X
  } else if (is.character(auxdesign)) {
    if (auxdesign == "fitted.values") {
      Z <- t(t(yhat))
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

  q <- ncol(Z) - 1
  n <- nrow(Z)
  auxresponse <- abs(e)
  auxres <- stats::lm.fit(Z, auxresponse)$residuals
  if (sigmaest == "main") {
    sigma_hatsq <- sum(e ^ 2) / n
  } else if (sigmaest == "auxiliary") {
    sigma_hatsq <- sum(auxres ^ 2) / n
  }

  teststat <- (sum(auxresponse ^ 2) - n * mean(auxresponse) ^ 2 -
                 sum(auxres ^ 2)) / (sigma_hatsq * (1 - 2 / pi))
  if (statonly) return(teststat)

  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
               null.value = "Homoskedasticity",
               alternative = "greater"), class = "htest")
  broom::tidy(rval)
}
