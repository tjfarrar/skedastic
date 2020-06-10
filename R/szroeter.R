#' Szroeter's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Szroeter78;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model.
#'
#' The test entails putting the data rows in increasing order of
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
#'    Default \code{NULL} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#'
#' @inheritParams breusch_pagan
#' @inheritParams carapeto_holt
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
#' szroeter(mtcars_lm, deflator = "qsec")
#'

szroeter <- function (mainlm, deflator = NULL, h = SKH,
                        qfmethod = "imhof") {

  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
    p <- ncol(X)
    hasintercept <- columnof1s(X)
    y <- stats::model.response(stats::model.frame(mainlm))
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
    p <- ncol(X)
    hasintercept <- columnof1s(X)
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)

  if (is.numeric(deflator) && deflator == as.integer(deflator)) {
    if (hasintercept[[1]] && deflator == 1) {
      stop("deflator cannot be the model intercept")
    } else if (deflator > p) {
      stop("`deflator` is not the index of a column of design matrix")
    }
  } else if (is.character(deflator)) {
    if (deflator == "(Intercept)") {
      stop("deflator cannot be the model intercept")
    } else if (!deflator %in% colnames(X)) {
      stop("`deflator` is not the name of a column of design matrix")
    }
  } else if (!is.null(deflator)) stop("`deflator` must be integer or character")

  if (!is.null(deflator)) {
    y <- y[order(X[, deflator])]
    X <- X[order(X[, deflator]), , drop = FALSE]
  }
  e <- stats::lm.fit(X, y)$residuals
  M <- fastM(X, n)
  Delta <- diag(h(1:n))

  teststat <- as.double((t(e) %*% Delta %*% e) / crossprod(e))
  pval <- pvalRQF(r = teststat, A = M %*% Delta %*% M,
                  B = M, algorithm = qfmethod, lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, p.value = pval,
            null.value = 1,
            alternative = "greater"), class = "htest")
  broom::tidy(rval)
}
