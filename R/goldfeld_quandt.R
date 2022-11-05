#' Goldfeld-Quandt Tests for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods (parametric and nonparametric) of
#'    \insertCite{Goldfeld65;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' @details The parametric test entails putting the data rows in increasing
#'    order of some specified deflator (one of the explanatory variables). A
#'    specified proportion of the most central observations (under this
#'    ordering) is removed, leaving a subset of lower observations and a subset
#'    of upper observations. Separate OLS regressions are fit to these two
#'    subsets of observations (using all variables from the original model).
#'    The test statistic is the ratio of the sum of squared residuals from the
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
#'    Default \code{NA} means the data will be left in its current order
#'    (e.g. in case the existing index is believed to be associated with
#'    error variance).
#' @param prop_central A double specifying the proportion of central
#'    observations to exclude from the F test (when \code{method} is
#'    \code{"parametric"} only). \code{\link[base]{round}} is
#'    used to ensure the number of central observations is an integer. The
#'    value must be small enough to allow the two auxiliary regressions to
#'    be fit; otherwise an error is thrown. Defaults to \code{1 / 3}.
#' @param group1prop A double specifying the proportion of remaining
#'    observations \emph{(after excluding central observations)} to allocate
#'    to the first group. The default value of \code{1 / 2} means that an
#'    equal number of observations is assigned to the first and second groups.
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
#'    \code{NA} (the default), probabilities are calculated within the
#'    function by calling \code{ppeak}. The user can improve computational
#'    performance of the test (for instance, when the test is being used
#'    repeatedly in a simulation) by pre-specifying the exact probability
#'    distribution of the number of peaks using this argument, e.g. by
#'    calling the \eqn{n}th element of \code{\link{dpeakdat}} (or \eqn{(n-p)}th
#'    element, if BLUS residuals are used).
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
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[lmtest:gqtest]{lmtest::gqtest}}, another implementation
#'    of the Goldfeld-Quandt Test (parametric method only).
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", prop_central = 0.25)
#' # This is equivalent to lmtest::gqtest(mtcars_lm, fraction = 0.25, order.by = mtcars$qsec)
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", method = "nonparametric",
#'  restype = "blus")
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", prop_central = 0.25, alternative = "two.sided")
#' goldfeld_quandt(mtcars_lm, deflator = "qsec", method = "nonparametric",
#'  restype = "blus", alternative = "two.sided")

goldfeld_quandt <- function(mainlm, method = c("parametric", "nonparametric"),
                    deflator = NA, prop_central = 1 / 3, group1prop = 1 / 2,
                    alternative = c("greater", "less", "two.sided"),
                    prob = NA, twosidedmethod = c("doubled", "kulinskaya"),
                    restype = c("ols", "blus"), statonly = FALSE, ...) {

  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  method <- match.arg(method, c("parametric", "nonparametric"))
  restype <- match.arg(restype, c("ols", "blus"))

  processmainlm(m = mainlm, needy = (method == "parametric"))

  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)

  checkdeflator(deflator, X, p, hasintercept[[1]])

  if (method == "parametric") {

    theind <- gqind(n, p, prop_central, group1prop)

    if (!is.na(deflator) && !is.null(deflator)) {
      if (!is.na(suppressWarnings(as.integer(deflator)))) {
        deflator <- as.integer(deflator)
      }
      y <- y[order(X[, deflator])]
      X <- X[order(X[, deflator]), , drop = FALSE]
    }
    thedf2 <- (length(theind[[2]]) - p)
    thedf1 <- (length(theind[[1]]) - p)
    S2sq <- sum(stats::lm.fit(X[theind[[2]], , drop = FALSE],
                              y[theind[[2]]])$residuals ^ 2) / thedf2
    S1sq <- sum(stats::lm.fit(X[theind[[1]], , drop = FALSE],
                              y[theind[[1]]])$residuals ^ 2) / thedf1
    teststat <- S2sq / S1sq
    if (statonly) return(teststat)
    names(thedf1) <- "df1"
    if (alternative == "greater") {
      pval <- stats::pf(teststat, df1 = thedf1, df2 = thedf2, lower.tail = FALSE)
    } else if (alternative == "less") {
      pval <- stats::pf(teststat, df1 = thedf1, df2 = thedf2, lower.tail = TRUE)
    } else if (alternative == "two.sided") {
      pval <- twosidedpval(q = teststat, CDF = stats::pf,
                       locpar = thedf2 / (thedf2 - 2),
                       method = twosidedmethod, continuous = TRUE,
                       df1 = thedf1, df2 = thedf2,
                       lower.tail = TRUE)
    }
    fullmethod <- "Goldfeld-Quandt F Test"
  } else if (method == "nonparametric") {
    if (restype == "ols") {
      absres <- abs(e)
      newn <- n
    } else if (restype == "blus") {
      absres <- abs(blus(mainlm = list("e" = e, "X" = X), ...))
      newn <- n - p
    }

    if (!is.na(deflator) && !is.null(deflator)) {
      if (!is.na(suppressWarnings(as.integer(deflator)))) {
        deflator <- as.integer(deflator)
      }
      absres <- absres[order(X[, deflator])]
    }

    teststat <- countpeaks(absres[!is.na(absres)])
    if (statonly) return(teststat)

    thedf1 <- NULL
    if (is.na(prob[1]) || is.null(prob)) {
      if (alternative == "greater") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = FALSE, usedata = (newn <= 1000))
      } else if (alternative == "less") {
        pval <- ppeak(k = teststat, n = newn, lower.tail = TRUE, usedata = (newn <= 1000))
      } else if (alternative == "two.sided") {
          peakmean <- sum(0:(newn - 1) * dpeak(k = 0:(newn - 1), n = newn,
                                               usedata = (newn <= 1000)))
          pval <- twosidedpval(q = teststat, locpar = peakmean, CDF = ppeak,
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
          pval <- twosidedpval(q = teststat, locpar = peakmean, CDF = pfunc,
                           method = twosidedmethod, continuous = FALSE)
      }
    }
    fullmethod <- "Goldfeld-Quandt Peaks Test"
  } else stop("Invalid `method` argument")

  rval <- structure(list(statistic = teststat, p.value = pval, parameter = thedf1,
               null.value = "Homoskedasticity", alternative = alternative,
               method = fullmethod), class = "htest")
  broom::tidy(rval)
}
