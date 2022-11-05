#' Szroeter's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Szroeter78;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model.
#'
#' @details The test entails putting the data rows in increasing order of
#'    some specified deflator (e.g., one of the explanatory variables) that
#'    is believed to be related to the error variance by some non-decreasing
#'    function. The test statistic is a ratio of quadratic forms in the OLS
#'    residuals. It is a right-tailed test.
#'
#' @param h A non-decreasing function taking as its argument the index
#'    \code{i} of observations from 1 to \eqn{n}. Defaults to \code{SKH},
#'    which is equivalent to \eqn{h(i)=2(1-\cos \frac{\pi i}{n+1})}.
#'    The function must be able to take a vector argument of length \code{n}.
#' @param deflator Either a character specifying a column name from the
#'    design matrix of \code{mainlm} or an integer giving the index of a
#'    column of the design matrix. This variable is suspected to be
#'    related to the error variance under the alternative hypothesis.
#'    \code{deflator} may not correspond to a column of 1's (intercept).
#'    Default \code{NA} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#'
#' @inheritParams breusch_pagan
#' @inheritParams carapeto_holt
#'
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' szroeter(mtcars_lm, deflator = "qsec")
#'

szroeter <- function (mainlm, deflator = NA, h = SKH,
                        qfmethod = "imhof", statonly = FALSE) {

  processmainlm(m = mainlm, needy = FALSE)

  n <- nrow(X)
  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  checkdeflator(deflator, X, p, hasintercept[[1]])

  if (!is.na(deflator) && !is.null(deflator)) {
    if (!is.na(suppressWarnings(as.integer(deflator)))) {
      deflator <- as.integer(deflator)
    }
    e <- e[order(X[, deflator])]
    X <- X[order(X[, deflator]), , drop = FALSE]
  }
  M <- fastM(X, n)
  Delta <- diag(h(1:n))

  teststat <- as.double((t(e) %*% Delta %*% e) / crossprod(e))
  if (statonly) return(teststat)
  pval <- pRQF(r = teststat, A = M %*% Delta %*% M,
                  B = M, algorithm = qfmethod, lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, p.value = pval,
            null.value = 1,
            alternative = "greater"), class = "htest")
  broom::tidy(rval)
}
