#' Glejser Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Glejser69;textual}{skedastic} for testing for "multiplicative"
#'    heteroskedasticity in a linear regression model.
#'    \insertCite{Mittelhammer00;textual}{skedastic} gives the
#'    formulation of the test used here.
#'
#' Glejser's Test entails fitting an auxiliary regression model in
#'    which the response variable is the absolute residual from the original
#'    model and the design matrix \eqn{Z} consists of one or more exogenous
#'    variables that are suspected of being related to the error variance.
#'    In the absence of prior information on a possible choice of \eqn{Z},
#'    one would typically use the explanatory variables from the original model.
#'    Under the null hypothesis of homoskedasticity, the distribution of the
#'    test statistic is asymptotically chi-squared with \code{parameter} degrees
#'    of freedom. The test is right-tailed.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso the description of the test in
#'   \href{http://www.econometrics.com/intro/testhet.htm}{SHAZAM} software
#'   (which produces identical results).
#'
#' @examples
#' n <- 20
#' p <- 4
#' set.seed(9586)
#' X <- matrix(data = runif(n * (p - 1)), nrow = n, ncol = p - 1)
#' # Response values generated under homoskedasticity
#' y_H0 <- rnorm(n, mean = 1 + rowSums(X), sd = 1)
#' glejser(lm(y_H0 ~ X))
#'# Response values generated under heteroskedasticity associated with X
#' y_HA <- rnorm(n, mean = 1 + rowSums(X), sd = rowSums(X ^ 2))
#' harvey(lm(y_HA ~ X))
#'

glejser <- function (mainlm, auxdesign = NULL) {

  if (class(mainlm) == "lm") {
    X <- model.matrix(mainlm)
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x), is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, ]
    }
    mainlm <- lm.fit(X, y)
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
  auxresponse <- abs(mainlm$residuals)
  auxres <- lm.fit(Z, auxresponse)$residuals
  sigma_hatsq <- sum(mainlm$residuals ^ 2) / n

  teststat <- (sum(auxresponse ^ 2) - n * mean(auxresponse) ^ 2 - sum(auxres ^ 2)) / (sigma_hatsq * (1 - 2 / pi))
  pval <- 1 - pchisq(teststat, df = p)

  rval <- structure(list(statistic = teststat, parameter = p, p.value = pval,
               null.value = "Homoskedasticity",
               alternative = "Heteroskedasticity"), class = "htest")
  broom::tidy(rval)
}
