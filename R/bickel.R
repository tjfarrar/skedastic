#' Bickel's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Bickel78;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model, with or without the
#'    scale-invariance modification of
#'    \insertCite{Carroll81;textual}{skedastic}.
#'
#' Bickel's Test is a robust extension of Anscombe's Test
#'    \insertCite{Anscombe61}{skedastic} in which the OLS residuals and
#'    estimated standard error are replaced with an \eqn{M} estimator. Under
#'    the null hypothesis of homoskedasticity, the distribution of the test
#'    statistic is asymptotically standard normally distributed. The test is
#'    two-tailed.
#' @param fitmethod A character indicating the method to be used to fit the
#'    regression model. This can be "OLS" for ordinary least squares (the
#'    default) or "robust" in which case a robust fitting method is called
#'    from \code{\link[MASS]{rlm}}.
#' @param a A character argument specifying the name of a function to be
#'    applied to the fitted values, or an integer \eqn{m} in which case the
#'    function applied is \eqn{f(x) = x^m}. Defaults to \code{"identity"} for
#'    \code{\link[base]{identity}}.
#' @param b A character argument specifying a function to be applied to the
#'    residuals. Defaults to Huber's function squared, as recommended by
#'    \insertCite{Carroll81;textual}{skedastic}. Currently the only supported
#'    functions are \code{"hubersq"} (for Huber's function squared) and
#'    \code{"tanhsq"} (for \eqn{b(x)=\mathrm{tanh}(x)^2}.)
#' @param scale_invariant A logical indicating whether the scale-invariance
#'    modification proposed by \insertCite{Carroll81;textual}{skedastic}
#'    should be implemented. Defaults to \code{TRUE}.
#' @param k A double argument specifying a parameter for Huber's function
#'    squared; used only if \code{b == "hubersq"}. This is not to be confused
#'    with the argument \code{k2} that could be passed to
#'    \code{\link[MASS]{rlm}} if the regression is fitted using robust methods.
#'    \code{k} defaults to 1.345.
#' @param ... Optional arguments to be passed to \code{\link[MASS]{rlm}}
#'
#' @inheritParams breusch_pagan
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso discussions of this test in
#'   \insertCite{Carroll81;textual}{skedastic} and
#'   \insertCite{Ali84;textual}{skedastic}.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' bickel(mtcars_lm)
#' bickel(mtcars_lm, fitmethod = "rlm")
#' bickel(mtcars_lm, scale_invariant = FALSE)
#'

bickel <- function (mainlm, fitmethod = c("lm", "rlm"),
                    a = "identity", b = c("hubersq", "tanhsq"),
                    scale_invariant = TRUE, k = 1.345, ...) {

  fitmethod <- match.arg(fitmethod, c("lm", "rlm"))
  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
    y <- stats::model.response(stats::model.frame(mainlm))
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x), is.nan(x),
                                                           is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, drop = FALSE]
    }
  }
  if (fitmethod == "lm") {
    mainlm <- stats::lm.fit(X, y)
  } else if (fitmethod == "rlm") {
    mainlm <- MASS::rlm(X, y)
  }

  p <- length(mainlm$coefficients)
  n <- nrow(X)
  yhat <- mainlm$fitted.values
  e <- mainlm$residuals
  afunc <- suppressWarnings(ifelse(is.na(as.double(a)), get(a),
                                   function(x) `^`(x, as.double(a))))
  b <- match.arg(b, c("hubersq", "tanhsq"))
  if (b == "hubersq") {
    bfunc <- function(x) {
      vapply(x, function(y) ifelse(abs(y) <= k, y ^ 2, k ^ 2), NA_real_)
    }
  } else if (b == "tanhsq") {
    bfunc <- function(x) tanh(x) ^ 2
  }

  if (scale_invariant) {
    sigma_hat <- stats::median(abs(e)) / stats::qnorm(0.75)
    e <- e / sigma_hat
  }

  sigma_b_sq <- 1 / (n - p) * sum((afunc(yhat) - mean(afunc(yhat))) ^ 2) *
    sum((bfunc(e) - mean(bfunc(e))) ^ 2)
  teststat <- sum(afunc(yhat) - mean(afunc(yhat)) * bfunc(e)) / sigma_b_sq
  pval <- 2 * stats::pnorm(abs(teststat), lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, p.value = pval,
               null.value = 0,
               alternative = "two.sided"), class = "htest")
  broom::tidy(rval)
}
