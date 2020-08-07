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
#' \eqn{n \times q} matrix consisting of \eqn{q} exogenous variables, appended
#' to a column of ones.
#' The test statistic is half the residual sum of squares from this auxiliary
#' regression. Under the null hypothesis of homoskedasticity, the distribution
#' of the test statistic is asymptotically chi-squared with \eqn{q} degrees of
#' freedom. The test is right-tailed.
#' @param hetfun A character describing the form of \eqn{w(\cdot)}, the error
#'    variance function under the heteroskedastic alternative. Possible values
#'    are \code{"mult"} (the default), corresponding to
#'    \eqn{w(Z_i,\lambda)=\exp\left\{\sum_{j=1}^{q}\lambda_j Z_{ij}\right\}},
#'    \code{"add"}, corresponding to
#'    \eqn{w(Z_i,\lambda)=\left(1+\sum_{j=1}^{q} \lambda_j Z_{ij}\right)^2}, and
#'    \code{"logmult"}, corresponding to
#'    \eqn{w(Z_i,\lambda)=\exp\left\{\sum_{j=1}^{q}\lambda_j \log Z_{ij}\right\}}.
#'    The multiplicative and log-multiplicative cases are considered in
#'    \insertCite{Cook83;textual}{skedastic}; the additive case is discussed,
#'    \emph{inter alia}, by \insertCite{Griffiths86;textual}{skedastic}.
#'    Results for the additive and multiplicative models are identical for this
#'    test. Partial matching is used.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[car:ncvTest]{car::ncvTest}}, which implements the same
#' test. Calling \code{car::ncvTest} with \code{var.formula} argument omitted
#' is equivalent to calling \code{skedastic::cook_weisberg} with
#' \code{auxdesign = "fitted.values", hetfun = "additive"}. Calling
#' \code{car::ncvTest} with \code{var.formula = ~ X} (where \code{X} is the
#' design matrix of the linear model with the intercept column omitted) is
#' equivalent to calling \code{skedastic::cook_weisberg} with default
#' \code{auxdesign} and \code{hetfun} values. The
#' \code{hetfun = "multiplicative"} option has no equivalent in
#' \code{car::ncvTest}.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' cook_weisberg(mtcars_lm)
#' cook_weisberg(mtcars_lm, auxdesign = "fitted.values", hetfun = "logmult")
#'

cook_weisberg <- function(mainlm, auxdesign = NA,
                  hetfun = c("mult", "add", "logmult"), statonly = FALSE) {

  hetfun <- match.arg(hetfun, c("mult", "add", "logmult"))

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
  if (hasintercept[[1]]) {
    Z <- Z[, -hasintercept[[2]], drop = FALSE]
  }

  q <- ncol(Z)
  n <- nrow(Z)

  if (hetfun == "mult") {
    Z <- cbind(1, Z)
  } else if (hetfun == "logmult") {
    Z <- cbind(1, log(Z))
  } else if (hetfun == "add") {
    Z <- cbind(1, 2 * Z)
  } else stop("Invalid hetfun argument")

  sigma_hatsq <- sum(e ^ 2) / n
  std_res_sq <- e ^ 2 / sigma_hatsq
  auxres <- stats::lm.fit(Z, std_res_sq)$residuals
  teststat <- (sum(std_res_sq ^ 2) - n * mean(std_res_sq) ^ 2 - sum(auxres ^ 2)) / 2
  if (statonly) return(teststat)
  method <- hetfun
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = method),
                         class = "htest")
  broom::tidy(rval)
}
