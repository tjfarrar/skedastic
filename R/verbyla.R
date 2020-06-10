#' Verbyla's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the residual maximum likelihood test of
#'    \insertCite{Verbyla93;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' Verbyla's Test entails fitting a generalised auxiliary regression model
#' in which the response variable is the vector of standardised squared
#' residuals \eqn{e_i^2/\hat{\sigma}^2} from the original OLS
#' model and the design matrix is some function of \eqn{Z}, an
#' \eqn{n \times q} matrix consisting of \eqn{q} exogenous variables, appended
#' to a column of ones.
#' The test statistic is half the residual sum of squares from this generalised
#' auxiliary regression. Under the null hypothesis of homoskedasticity, the
#' distribution of the test statistic is asymptotically chi-squared with
#' \eqn{q} degrees of freedom. The test is right-tailed.
#' @inheritParams breusch_pagan
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
#' verbyla(mtcars_lm)
#' verbyla(mtcars_lm, auxdesign = "fitted.values")
#'

verbyla <- function (mainlm, auxdesign = NULL) {

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
  } else if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x),
                                                    is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, drop = FALSE]
    }
    mainlm <- stats::lm.fit(X, y)
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
    Z <- Z[, -hasintercept[[2]], drop = FALSE]
  }

  p <- ncol(X)
  q <- ncol(Z)
  n <- nrow(Z)

  Z <- cbind(1, Z)

  M <- fastM(X, n)
  sigma_hatbar <- sum(mainlm$residuals ^ 2) / (n - p)
  term1 <- t(t(mainlm$residuals ^ 2 / sigma_hatbar - diag(M)))

  teststat <- as.double(1 / 2 * t(term1) %*% Z %*% solve(t(Z) %*% (M ^ 2) %*% Z)
                        %*% t(Z) %*% term1)

  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = "Verbyla"),
                         class = "htest")
  broom::tidy(rval)
}
