#' Estimate Covariance Matrix of Ordinary Least Squares Estimators
#'    Using Error Variance Estimates from an Auxiliary Variance Model
#'
#' The function simply calculates
#'    \deqn{\mathrm{Cov}{\hat{\beta}}=(X'X)^{-1}X'\hat{\Omega}X(X'X)^{-1}},
#'    where \eqn{X} is the design matrix of a linear regression model and
#'    \eqn{\hat{\Omega}} is an estimate of the diagonal variance-covariance
#'    matrix of the random errors, whose diagonal elements have been
#'    obtained from an auxiliary variance model fit with \code{alvm.fit}
#'    or \code{anlvm.fit}.
#'
#' @param as_matrix A logical. If \code{TRUE} (the default), a
#'   \eqn{p \times p} matrix is returned, where \eqn{p} is the
#'   number of columns in \eqn{X}. Otherwise, a numeric vector of length
#'   \eqn{p} is returned.
#'
#' @inheritParams avm.fwls
#'
#' @return Either a numeric matrix or a numeric vector, whose (diagonal)
#'   elements are \eqn{\widehat{\mathrm{Var}}(\hat{\beta}_j)},
#'   \eqn{j=1,2,\ldots,p}.
#'
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link{alvm.fit}}, \code{\link{anlvm.fit}},
#'    \code{\link{avm.fwls}}. If a matrix is returned, it can be
#'    passed to \code{\link[lmtest]{coeftest}} for implementation
#'    of a quasi-\eqn{t}-test of significance of the \eqn{\beta} coefficients.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' myalvm <- alvm.fit(mainlm = mtcars_lm, model = "linear",
#'    varselect = "qgcv.linear")
#' myvcov <- avm.vcov(myalvm)
#' lmtest::coeftest(mtcars_lm, vcov. = myvcov)
#'

avm.vcov <- function(object, as_matrix = TRUE) {

  X <- stats::model.matrix(object$ols)

  var.beta.hats <- Rfast::spdinv(crossprod(X)) %*% t(X) %*%
      diag(object$var.est) %*% X %*% Rfast::spdinv(crossprod(X))
  if (!as_matrix) {
    diag(var.beta.hats)
  } else {
    var.beta.hats
  }
}
