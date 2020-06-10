#' Honda's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two-sided LM Test of
#'    \insertCite{Honda89;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The test assumes that heteroskedasticity, if present, would be either of
#'    the additive form \eqn{\sigma_i^2 = \sigma^2(1+\theta z_i)} or of the
#'    multiplicative form \eqn{\sigma_i^2 = \sigma^2 e^{\theta z_i}}, where
#'    where \eqn{z_i} is a deflator (a nonstochastic variable
#'    suspected of being related to the error variance), \eqn{\sigma^2} is
#'    some unknown constant, and \eqn{\theta} is an unknown parameter
#'    representing the degree of heteroskedasticity. Since the test
#'    statistic \eqn{Q=\frac{e' A_0 e}{e'e}} is a ratio of quadratic forms
#'    in the errors, the Imhof algorithm is used to compute \eqn{p}-values.
#'    Since the null distribution is in general asymmetrical, the two-sided
#'    \eqn{p}-value is computed using the conditional method of
#'    \insertCite{Kulinskaya08;textual}{skedastic}, setting \eqn{A=1}.
#'
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
#' honda(mtcars_lm, deflator = "qsec")
#'

honda <- function (mainlm, deflator = NULL, alternative = c("two.sided",
                    "greater", "less"), twosidedmethod =
                     c("doubled", "kulinskaya"), qfmethod = "imhof") {

  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))
  if (class(mainlm) == "lm") {
    X <- stats::model.matrix(mainlm)
    hasintercept <- columnof1s(X)
    e <- mainlm$residuals
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
    hasintercept <- columnof1s(X)
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
    e <- stats::lm.fit(X, y)$residuals
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

  M <- fastM(X, n)
  if (is.null(deflator)) {
    A0 <- diag(1:n)
  } else {
    A0 <- diag(X[, deflator])
  }

  teststat <- as.double((t(e) %*% A0 %*% e) / crossprod(e))
  if (alternative == "greater") {
    pval <- pvalRQF(r = teststat, A = M %*% A0 %*% M,
                    B = M, algorithm = qfmethod, lower.tail = FALSE)
  } else if (alternative == "less") {
    pval <- pvalRQF(r = teststat, A = M %*% A0 %*% M,
                    B = M, algorithm = qfmethod, lower.tail = TRUE)
  } else if (alternative == "two.sided") {
    pval <- twosidedpval(q = teststat, CDF = pvalRQF,
                         method = twosidedmethod, Aloc = 1, continuous = TRUE,
                         A = M %*% A0 %*% M, B = M, algorithm = qfmethod,
                         lower.tail = TRUE)
  }

  rval <- structure(list(statistic = teststat, p.value = pval,
            null.value = 1,
            alternative = alternative), class = "htest")
  broom::tidy(rval)
}
