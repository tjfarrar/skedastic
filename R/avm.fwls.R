#' Apply Feasible Weighted Least Squares to a Linear Regression Model
#'
#' This function applies feasible weighted least squares (FWLS) to a
#'    linear regression model using error variance estimates obtained
#'    from an auxiliary linear variance model fit using \code{alvm.fit}
#'    or from an auxiliary nonlinear variance model fit using
#'    \code{anlvm.fit}.
#'
#' The function simply calculates
#'    \deqn{\hat{\beta}=(X'\hat{\Omega}^{-1}X)^{-1}X'\hat{\Omega}^{-1}y},
#'    where \eqn{X} is the design matrix, \eqn{y} is the response vector, and
#'    \eqn{\hat{\Omega}} is the diagonal variance-covariance matrix of the
#'    random errors, whose diagonal elements have been estimated by an
#'    auxiliary variance model.
#'
#' @param object Either an object of class \code{"alvm.fit"} or an object
#'    of class \code{"anlvm.fit"}
#' @param fastfit A logical. If \code{FALSE} (the default), the linear
#'    regression model is fit using \code{\link[stats]{lm}}; otherwise,
#'    using \code{\link[stats]{lm.wfit}}
#'
#' @return Either an object of \code{\link[base]{class}} \code{"lm"}
#'    (if \code{fastfit} is \code{FALSE}) or otherwise a generic
#'    list object
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link{alvm.fit}}, \code{\link{anlvm.fit}},
#'    \code{\link{avm.vcov}}
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' myalvm <- alvm.fit(mainlm = mtcars_lm, model = "linear",
#'    varselect = "qgcv.linear")
#' myfwls <- avm.fwls(myalvm)
#' cbind(coef(mtcars_lm), coef(myfwls))
#'

avm.fwls <- function(object, fastfit = FALSE) {
  X <- stats::model.matrix(object$ols)
  y <- stats::model.response(stats::model.frame(object$ols))

  if (fastfit) {
    stats::lm.wfit(x = X, y = y, w = 1 / object$var.est)
  } else {
    stats::lm(y ~ 0 + X, weights = 1 / object$var.est)
  }
}
