#' Carapeto-Holt Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods (parametric and nonparametric) of
#'    \insertCite{Carapeto03;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The test is based on the methodology of
#'    \insertCite{Goldfeld65;textual}{skedastic} but does not require any
#'    auxiliary regression. It entails ordering the observations by some
#'    suspected deflator (one of the explanatory variables) in such a way
#'    that, under the alternative hypothesis, the observations would now
#'    be arranged in decreasing order of error variance. A specified proportion
#'    of the most central observations (under this ordering) is removed,
#'    leaving a subset of lower observations and a subset of upper
#'    observations. The test statistic is then computed as a ratio of quadratic
#'    forms corresponding to the sums of squared residuals of the upper and
#'    lower observations respectively. \eqn{p}-values are computed by the
#'    Imhof algorithm in \code{\link{pvalRQF}}.
#'
#' @param deflator Either a character specifying a column name from the
#'    design matrix of \code{mainlm} or an integer giving the index of a
#'    column of the design matrix. This variable is suspected to be
#'    related to the error variance under the alternative hypothesis.
#'    \code{deflator} may not correspond to a column of 1's (intercept).
#'    Default \code{NULL} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#' @param prop_central A double specifying the proportion of central
#'    observations to exclude when comparing the two subsets of observations.
#'    \code{\link[base]{round}} is used to ensure the number of central
#'    observations is an integer. Defaults to \eqn{\frac{1}{3}}.
#' @param qfmethod A character, either \code{"imhof"}, \code{"davies"}, or
#'    \code{"integrate"}, corresponding to the \code{algorithm} argument
#'    of \code{\link[skedastic]{pvalRQF}}. The default is \code{"imhof"}.
#' @param alternative A character specifying the form of alternative
#'    hypothesis. If it is suspected that the error variance is positively
#'    associated with the deflator variable, \code{"greater"}. If it is
#'    suspected that the error variance is negatively associated with deflator
#'    variable, \code{"less"}. If no information is available on the suspected
#'    direction of the association, \code{"two.sided"}. Defaults to
#'    \code{"greater"}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams goldfeld_quandt
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
#' carapeto_holt(mtcars_lm, deflator = "qsec", prop_central = 0.25)
#' # Same as previous example
#' mtcars_list <- list("y" = mtcars$mpg, "X" = cbind(1, mtcars$wt, mtcars$qsec, mtcars$am))
#' carapeto_holt(mtcars_list, deflator = 3, prop_central = 0.25)
#'

carapeto_holt <- function(mainlm, deflator = NULL, prop_central = 1 / 3,
                  group1prop = 1 / 2, qfmethod = "imhof",
                  alternative = c("greater", "less", "two.sided"),
                  twosidedmethod = c("doubled", "kulinskaya"),
                  statonly = FALSE) {

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))

  processmainlm(m = mainlm, needy = FALSE)

  hasintercept <- columnof1s(X)
  if (class(mainlm) == "list") {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column
                                       of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)

  checkdeflator(deflator, X, p, hasintercept[[1]])

  theind <- gqind(n, prop_central, group1prop)

  if (!is.null(deflator)) {
    e <- e[order(X[, deflator], decreasing = TRUE)]
    X <- X[order(X[, deflator], decreasing = TRUE), , drop = FALSE]
  }

  M <- fastM(X, n)
  Istar_hi <- diag(x = 0, nrow = n)
  Istar_lo <- Istar_hi
  diag(Istar_hi)[theind[[2]]] <- 1
  diag(Istar_lo)[theind[[1]]] <- 1

  if (hasintercept[[1]]) {
    teststat <- as.double((t(e) %*% Istar_hi %*% e) /
                            (t(e) %*% Istar_lo %*% e))
    if (statonly) return(teststat)

    if (alternative == "greater") {
      pval <- pvalRQF(r = teststat, A = M %*% Istar_hi %*% M,
                      B = M %*% Istar_lo %*% M, algorithm = qfmethod,
                      lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- pvalRQF(r = teststat, A = M %*% Istar_hi %*% M,
                      B = M %*% Istar_lo %*% M, algorithm = qfmethod,
                      lower.tail = TRUE)
    } else if (alternative == "two.sided") {
      pval <- twosidedpval(q = teststat, CDF = pvalRQF, method = twosidedmethod,
              continuous = TRUE, Aloc = 1, A = M %*% Istar_hi %*% M,
              B = M %*% Istar_lo %*% M, algorithm = qfmethod, lower.tail = TRUE)
    }
  } else {
    A <- diag(n) - matrix(data = 1 / n, nrow = n, ncol = n)
    teststat <- as.double((t(A %*% e) %*% Istar_hi %*% A %*% e) /
      (t(A %*% e) %*% Istar_lo %*% A %*% e))
    if (statonly) return(teststat)

    if (alternative == "greater") {
      pval <- pvalRQF(r = teststat, A = M %*% A %*% Istar_hi %*% A %*% M,
                      B = M %*% A %*% Istar_lo %*% A %*% M, algorithm = qfmethod,
                      lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- pvalRQF(r = teststat, A = M %*% A %*% Istar_hi %*% A %*% M,
                      B = M %*% A %*% Istar_lo %*% A %*% M, algorithm = qfmethod,
                      lower.tail = TRUE)
    } else if (alternative == "two.sided") {
      pval <- twosidedpval(q = teststat, CDF = pvalRQF,
                      method = twosidedmethod, Aloc = 1, continuous = TRUE,
                       A = M %*% A %*% Istar_hi %*% A %*% M,
                       B = M %*% A %*% Istar_lo %*% A %*% M,
                       algorithm = qfmethod, lower.tail = TRUE)
    }
  }

  rval <- structure(list(statistic = teststat, p.value = pval, parameter = NULL,
               null.value = 1, alternative = alternative),
               class = "htest")
  broom::tidy(rval)
}
