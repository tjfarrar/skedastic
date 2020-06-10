#' Evans-King Tests for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods of
#'    \insertCite{Evans88;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The test entails putting the data rows in increasing order of
#'    some specified deflator (e.g., one of the explanatory variables) that
#'    is believed to be related to the error variance by some non-decreasing
#'    function. There are two statistics that can be used, corresponding to
#'    the two values of the \code{method} argument. In both cases the test
#'    statistic can be expressed as a ratio of quadratic forms in the errors,
#'    and thus the Imhof algorithm is used to compute \eqn{p}-values. Both
#'    methods involve a left-tailed test.
#'
#' @param method A character indicating which of the two tests derived in
#'    \insertCite{Evans88;textual}{skedastic} should be implemented.
#'    Possible values are \code{"GLS"} and \code{"LM"}; partial matching is
#'    used (which is not case-sensitive).
#' @param deflator Either a character specifying a column name from the
#'    design matrix of \code{mainlm} or an integer giving the index of a
#'    column of the design matrix. This variable is suspected to be
#'    related to the error variance under the alternative hypothesis.
#'    \code{deflator} may not correspond to a column of 1's (intercept).
#'    Default \code{NULL} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#' @param lambda_star A double; coefficient representing the degree of
#'    heteroskedasticity under the alternative hypothesis.
#'    \insertCite{Evans85;textual}{skedastic} suggests 2.5, 5, 7.5, and 10 as
#'    values to consider, and \insertCite{Evans88;textual}{skedastic} finds
#'    that 2.5 and 5 perform best empirically. This parameter is used only for
#'    the \code{"GLS"} method; the \code{"LM"} method represents the limiting
#'    case as \eqn{\lambda^\star \to 0}. Defaults to \code{5}.
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
#' @seealso \insertCite{Evans85;textual}{skedastic}, which already anticipates
#'    one of the tests.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' evans_king(mtcars_lm, deflator = "qsec", method = "GLS")
#' evans_king(mtcars_lm, deflator = "qsec", method = "LM")
#'

evans_king <- function (mainlm, method = c("GLS", "LM"), deflator = NULL,
                        lambda_star = 5, qfmethod = "imhof") {

  method <- match.arg(toupper(method), c("GLS", "LM"))

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

  if (method == "GLS") {
    tau <- ((1:n) - 1) / (n - 1)
    w <- (1 + lambda_star * tau) ^ (1 / 2)
    R <- diag(1 / w)
    Xstar <- t(sapply(1:n, function(i) X[i, ] / w[i], USE.NAMES = FALSE))
    Mstar <- fastM(Xstar, n)
    teststat <- as.double((t(e) %*% R %*% Mstar %*% R %*% e) / crossprod(e))
    pval <- pvalRQF(r = teststat, A = M %*% R %*% Mstar %*% R %*% M,
                    B = M, algorithm = qfmethod, lower.tail = TRUE)
  } else if (method == "LM") {
    lambda_star <- NULL
    N <- diag((n - 1:n) / (n - 1))
    teststat <- as.double((t(e) %*% N %*% e) / crossprod(e))
    pval <- pvalRQF(r = teststat, A = M %*% N %*% M,
                    B = M, algorithm = qfmethod, lower.tail = TRUE)
  } else stop("Invalid `method` argument")

  rval <- structure(list(statistic = teststat, p.value = pval,
            parameter = lambda_star, null.value = 1,
            alternative = "less", method = method), class = "htest")
  broom::tidy(rval)
}
