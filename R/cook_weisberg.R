#' Cook-Weisberg Score Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the score test of
#'    \insertCite{Cook83;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The Cook-Weisberg Score Test entails fitting an auxiliary regression model
#' in which the response variable is the vector of standardised squared
#' residuals \eqn{e_i^2/\hat{\sigma}^2} from the original OLS
#' model and the design matrix is some function of \eqn{Z}, an
#' \eqn{n \times q} matrix consisting of \eqn{q} exogenous variables.
#' The test statistic is half the residual sum of squares from this auxiliary
#' regression. Under the null hypothesis of homoskedasticity, the distribution
#' of the test statistic is asymptotically chi-squared with \eqn{q} degrees of
#' freedom. The test is right-tailed.
#' @param errorfun A character describing the functional form of
#'    the error variance under the heteroskedastic alternative. Possible values
#'    are \code{"additive"} (the default) and
#'    \code{"multiplicative"}, corresponding to the two cases considered in
#'    \insertCite{Cook83;textual}{skedastic}, or the name of a function in the
#'    environment (passed as a character). If the name of a function, it
#'    will be applied to auxiliary design element \eqn{z_{ij}} to obtain the
#'    corresponding element of \eqn{D}, according to the notation used in
#'    \insertCite{Cook83;textual}{skedastic}. The value \code{"additive"}
#'    corresponds to the function \code{identity} and \code{"multiplicative"}
#'    to the function \code{log}. Partial matching is NOT used.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[car:ncvTest]{car::ncvTest}}, which implements the same
#' test. Calling \code{car::ncvTest} with \code{var.formula} argument omitted
#' is equivalent to calling \code{skedastic::cook_weisberg} with
#' \code{auxdesign = "fitted.values", errorfun = "additive"}. Calling
#' \code{car::ncvTest} with \code{var.formula = ~ X} (where \code{X} is the
#' design matrix of the linear model with the intercept column omitted) is
#' equivalent to calling \code{skedastic::cook_weisberg} with default
#' \code{auxdesign} and \code{errorfun} values. The
#' \code{errorfun = "multiplicative"} option has no equivalent in
#' \code{car::ncvTest}.
#'
#' @examples
#' n <- 20
#' p <- 4
#' set.seed(9586)
#' X <- matrix(data = runif(n * (p - 1)), nrow = n, ncol = p - 1)
#' # Response values generated under homoskedasticity
#' y_H0 <- rnorm(n, mean = 1 + rowSums(X), sd = 1)
#' cook_weisberg(lm(y_H0 ~ X))
#' cook_weisberg(lm(y_H0 ~ X), auxdesign = "fitted.values")
#' cook_weisberg(lm(y_H0 ~ X), errorfun = "multiplicative")
#'# Response values generated under heteroskedasticity associated with X
#' y_HA <- rnorm(n, mean = 1 + rowSums(X), sd = rowSums(X ^ 2))
#' cook_weisberg(lm(y_HA ~ X))
#' cook_weisberg(lm(y_HA ~ X), auxdesign = "fitted.values")
#' cook_weisberg(lm(y_HA ~ X), errorfun = "multiplicative")
#'

cook_weisberg <- function (mainlm, auxdesign = NULL, errorfun = "additive") {

  if (class(mainlm) == "lm") {
    X <- model.matrix(mainlm)
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
  if (hasintercept[[1]]) {
    Z <- Z[, -hasintercept[[2]]]
  }

  q <- ncol(Z)
  n <- nrow(Z)

  if (errorfun == "additive") {
    Z <- cbind(1, Z)
  } else if (errorfun == "multiplicative") {
    Z <- cbind(1, log(Z))
  } else {
    Z <- cbind(1, get(errorfun)(Z))
  }

  sigma_hatsq <- sum(mainlm$residuals ^ 2) / n
  std_res_sq <- mainlm$residuals ^ 2 / sigma_hatsq
  auxres <- lm.fit(Z, std_res_sq)$residuals
  teststat <- (sum(std_res_sq ^ 2) - n * mean(std_res_sq) ^ 2 - sum(auxres ^ 2)) / 2
  method <- errorfun
  pval <- 1 - pchisq(teststat, df = q)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "Heteroskedasticity", method = method),
                         class = "htest")
  broom::tidy(rval)
}
