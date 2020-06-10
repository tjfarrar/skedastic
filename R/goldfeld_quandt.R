#' Goldfeld-Quandt Tests for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods (parametric and nonparametric) of
#'    \insertCite{Goldfeld65;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The parametric test entails putting the data rows in increasing order of
#'    some specified deflator (one of the explanatory variables). A specified
#'    proportion of the most central observations (under this ordering) is
#'    removed, leaving a subset of lower observations and a subset of upper
#'    observations. Separate OLS regressions are fit to these two subsets of
#'    observations (using all variables from the original model). The test
#'    statistic is the ratio of the sum of squared residuals from the
#'    'upper' model to the sum of squared residuals from the 'lower' model.
#'    Under the null hypothesis, the test statistic is exactly F-distributed
#'    with numerator and denominator degrees of freedom equal to
#'    \eqn{(n-c)/2 - p} where \eqn{n} is the number of observations in the
#'    original regression model, \eqn{c} is the number of central observations
#'    removed, and \eqn{p} is the number of columns in the design matrix (number of
#'    parameters to be estimated, including intercept).
#'
#' The nonparametric test entails putting the residuals of the linear model in
#'    increasing order of some specified deflator (one of the explanatory
#'    variables). The test statistic is the number of peaks, with the \eqn{j}th
#'    absolute residual \eqn{|e_j|} defined as a peak if \eqn{|e_j|\ge|e_i|}
#'    for all \eqn{i<j}. The first observation does not constitute a peak. If
#'    the number of peaks is large relative to the distribution of peaks under
#'    the null hypothesis, this constitutes evidence for heteroskedasticity.
#'
#' @param method A character indicating which of the two tests derived in
#'    \insertCite{Goldfeld65;textual}{skedastic} should be implemented.
#'    Possible values are "parametric" and "nonparametric". Default is
#'    "parametric". It is acceptable to specify only the first letter.
#' @param deflator Either a character specifying a column name from the
#'    design matrix of \code{mainlm} or an integer giving the index of a
#'    column of the design matrix. This variable is suspected to be
#'    related to the error variance under the alternative hypothesis.
#'    \code{deflator} may not correspond to a column of 1's (intercept).
#'    Default \code{NULL} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#' @param prop_central A double specifying the proportion of central
#'    observations to exclude from the F test (when \code{method} is
#'    \code{"parametric"} only). \code{\link[base]{round}} is
#'    used to ensure the number of central observations is an integer. The
#'    value must be small enough to allow the two auxiliary regressions to
#'    be fit; otherwise an error is thrown. Defaults to \eqn{\frac{1}{3}}.
#' @param alternative A character specifying the form of alternative
#'    hypothesis. If it is suspected that the
#'    error variance is positively associated with the deflator variable,
#'    \code{"greater"}. If it is suspected that the error variance is
#'    negatively associated with deflator variable, \code{"less"}. If no
#'    information is available on the suspected direction of the association,
#'    \code{"two.sided"}. Defaults to \code{"greater"}.
#' @param prob A vector of probabilities corresponding to values of the test
#'    statistic (number of peaks) from 0 to \eqn{n-1} inclusive (used
#'    only when \code{method} is \code{"nonparametric"}). If
#'    \code{NULL} (the default), probabilities are calculated within the
#'    function by calling \code{ppeak}. The user can improve computational
#'    performance of the test (for instance, when the test is being used
#'    repeatedly in a simulation) by pre-specifying the exact probability
#'    distribution of the number of peaks using this argument.
#' @param twosidedmethod A character indicating the method to be used to compute
#'    two-sided \eqn{p}-values for the parametric test when \code{alternative}
#'    is \code{"two.sided"}. The argument is passed to
#'    \code{\link{twosidedpval}} as its \code{method} argument.
#' @param restype A character specifying which residuals to use: \code{"ols"}
#'    for OLS residuals (the default) or the \code{"blus"} for
#'    \link[=blus]{BLUS} residuals. The advantage of using BLUS residuals is
#'    that, under the null hypothesis, the assumption that the random series
#'    is independent and identically distributed is met (whereas with OLS
#'    residuals it is not). The disadvantage of using BLUS residuals is that
#'    only \eqn{n-p} residuals are used rather than the full \eqn{n}. This
#'    argument is ignored if \code{method} is \code{"parametric"}.
#' @param ... Optional further arguments to pass to \code{\link{blus}}.
#'
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[lmtest:gqtest]{lmtest::gqtest}}, which performs the
#'    parametric version of the Goldfeld-Quandt Test. The `point` argument in
#'    that function allows the splitting point of data set into subsets to be
#'    other than the middle observation and thus allows the subsets to be of
#'    different sizes.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", prop_central = 0.25)
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", method = "nonparametric",
#'  restype = "blus")
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", prop_central = 0.25, alternative = "two.sided")
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", method = "nonparametric",
#'  restype = "blus", alternative = "two.sided")

goldfeld_quandt <- function (mainlm, method = "parametric", deflator = NULL,
                    prop_central = 1 / 3, alternative = c("greater", "less",
                    "two.sided"), prob = NULL, twosidedmethod = c("doubled",
                    "kulinskaya"), restype = c("ols", "blus"), ...) {

  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  method <- match.arg(method, c("parametric", "nonparametric"))
  restype <- match.arg(restype, c("ols", "blus"))

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
    if (method == "nonparametric") mainlm <- stats::lm.fit(X, y)
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

  if (method == "parametric") {

    k <- num_to_remove(n, prop_central)

    if (!is.null(deflator)) {
      y <- y[order(X[, deflator])]
      X <- X[order(X[, deflator]), , drop = FALSE]
    }
    ind_lo <- 1:((n - k) / 2)
    ind_hi <- ((n + k) / 2 + 1):n
    S_hi <- sum(stats::lm.fit(X[ind_hi, ], y[ind_hi])$residuals ^ 2)
    S_lo <- sum(stats::lm.fit(X[ind_lo, ], y[ind_lo])$residuals ^ 2)
    teststat <- S_hi / S_lo
    param <- (n - k) / 2 - p
    names(param) <- "df"
    if (alternative == "greater") {
      pval <- stats::pf(teststat, df1 = param, df2 = param, lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- stats::pf(teststat, df1 = param, df2 = param, lower.tail = TRUE)
    } else if (alternative == "two.sided") {
      pval <- twosidedpval(q = teststat, Aloc = 1, CDF = pf,
                       method = twosidedmethod, df1 = param, df2 = param)
    }
    fullmethod <- "Goldfeld-Quandt F Test"
  } else if (method == "nonparametric") {
    if (restype == "ols") {
      absres <- abs(mainlm$residuals)
      newn <- n
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm, ...))
      newn <- n - p
    }

    if (!is.null(deflator)) absres <- absres[order(X[, deflator])]

    teststat <- countpeaks(absres[!is.na(absres)])
    param <- NULL
    if (is.null(prob)) {
      if (alternative == "greater") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = FALSE, usedata = (newn <= 1000))
      } else if (alternative == "less") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = TRUE, usedata = (newn <= 1000))
      } else if (alternative == "two.sided") {
          peakmean <- sum(0:(newn - 1) * dpeak(k = 0:(newn - 1), n = newn,
                                               usedata = (newn <= 1000)))
          pval <- twosidedpval(q = teststat, Aloc = peakmean, CDF = ppeak,
                           method = twosidedmethod, continuous = FALSE,
                           n = newn, lower.tail = TRUE, usedata = (newn <= 1000))
      }
    } else {
      if (length(prob) != newn) stop("prob must be a vector of length equal to number of observations in series")
      if (alternative == "greater") {
        pval <- sum(prob[(teststat + 1):length(prob)])
      } else if (alternative == "less") {
        pval <- sum(prob[1:(teststat + 1)])
      } else if (alternative == "two.sided") {
          pfunc <- function(k) {
            sum(prob[1:(k + 1)])
          }
          peakmean <- sum(0:(newn - 1) * prob)
          pval <- twosidedpval(q = teststat, Aloc = peakmean, CDF = pfunc,
                           method = twosidedmethod, continuous = FALSE)
      }
    }
    fullmethod <- "Goldfeld-Quandt Peaks Test"
  } else stop("Invalid `method` argument")

  rval <- structure(list(statistic = teststat, p.value = pval, parameter = param,
               null.value = "Homoskedasticity", alternative = alternative,
               method = fullmethod), class = "htest")
  broom::tidy(rval)
}
