#' Dufour et al.'s Monte Carlo Test for Heteroskedasticity in a Linear
#' Regression Model
#'
#' This function implements the method of
#'    \insertCite{Dufour04;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model.
#'
#' The test implements a Monte Carlo procedure as follows. (1) The observed
#'    value of the test statistic \eqn{T_0} is computed using \code{hettest}.
#'    (2) \code{R} replications of the random error vector are generated from
#'    the distribution specified using \code{errorgen}. (3) \code{R}
#'    replications of the test statistic, \eqn{T_1,T_2,\ldots,T_R}, are
#'    computed from the generated error vectors. (4) The empirical
#'    \code{p}-value is computed as \eqn{\frac{\hat{G}_R(T_0)+1}{R+1}}, where
#'    \eqn{\hat{G}_R(x)=\sum_{j=1}^{R} 1_{T_j \ge x}}, \eqn{1_{\bullet}}
#'    being the indicator function. The test is right-tailed, regardless of the
#'    tailedness of \code{hettest}. Note that the heteroskedasticity
#'    test implemented by \code{hettest} must have a test statistic that is
#'    continuous and that is invariant with respect to nuisance parameters
#'    (\eqn{\sigma^2} and \eqn{\beta}). Note further that if \code{hettest}
#'    is \code{goldfeld_quandt} with \code{method} argument
#'    \code{"parametric"}, the replicated Goldfeld-Quandt \eqn{F} statistics
#'    are computed directly within this function rather than by calling
#'    \code{goldfeld_quandt}, due to some idiosyncratic features of this test.
#'    Note that, if \code{alternative} is set to \code{"two.sided"}, the
#'    one-sided \eqn{p}-value is doubled (\code{\link{twosidedpval}} cannot
#'    be used in this case).
#'
#' @param R An integer specifying the number of Monte Carlo replicates to
#'    generate. Defaults to \code{1000}.
#' @param hettest A function, or a character specifying the name of a function,
#'    that implements a heteroskedasticity test on a linear regression model.
#'    The function is called with the \code{statonly} argument set to
#'    \code{TRUE} to improve computational efficiency.
#' @param alternative The tailedness of the test whose statistic is computed by
#'    \code{hettest}; one of \code{"greater"} (the default), \code{"less"}, or
#'    \code{"two.sided"}.
#' @param errorgen A function, or a character specifying the name of a
#'    function, from which the random errors are to be generated. The function
#'    should correspond to a continuous probability distribution that has (or
#'    at least can have) a mean of 0. Defaults to \code{\link[stats]{rnorm}}.
#' @param errorparam An optional list of parameters to pass to \code{errorgen}.
#'    This argument is ignored if \code{errorgen} is
#'    \code{\link[stats]{rnorm}}, since \code{mean} must be set to 0, and
#'    \code{sd} is set to 1 because the heteroskedasticity test implemented by
#'    \code{hettest} is assumed to be scale invariant. If \code{errorgen} is
#'    not \code{rnorm}, \code{errorparam} should be chosen in such a way that
#'    the error distribution has a mean of 0.
#' @param seed An integer specifying a seed to pass to
#'    \code{\link[base]{set.seed}} for random number generation. This allows
#'    reproducibility of Monte Carlo results. A value of \code{NULL}
#'    results in not setting a seed.
#' @param ... Additional arguments to pass to \code{hettest}
#'
#' @inheritParams breusch_pagan
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object
#'    is not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' dufour_etal(mtcars_lm, hettest = breusch_pagan)
#'

dufour_etal <- function(mainlm, hettest, R = 1000L, alternative = c("greater",
                        "less", "two.sided"), errorgen = stats::rnorm,
                        errorparam = list(), seed = 1234, ...) {

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  if (is.character(hettest)) hettest <- get(hettest)
  if (is.character(errorgen)) errorgen <- get(errorgen)
  arguments <- list(...)
  invisible(list2env(arguments, envir = environment()))
  if ("alternative" %in% names(formals(hettest))) arguments$alternative <- alternative
  hettestchar <- deparse(substitute(hettest))
  processmainlm(m = mainlm)
  n <- length(e)

  if (hettestchar == "horn") {
    if (exists("exact", where = environment(), inherits = FALSE)) {
      if (exact) stop("This method is only available for tests with a continuous test statistic")
    } else {
      if (n - p <= 11) stop("This method is only available for tests with a continuous test statistic")
    }
  } else if (hettestchar %in% c("anscombe", "bickel")) {
    stop("This method is only available for tests that are invariant with respect
         to nuisance parameters and for which the test statistic can be computed
         based only on the design matrix and the OLS residuals.")
  }

  if (exists("deflator", where = environment(), inherits = FALSE)) {
    if (is.character(deflator))
      arguments$deflator <- which(colnames(X) == deflator)
  }

  statobs <- do.call(what = hettest,
              args = append(list("mainlm" = list("y" = y, "X" = X),
                                 "statonly" = TRUE), arguments))

  if (!(hettestchar == "goldfeld_quandt")) {
    M <- fastM(X, n)
    if (!is.null(seed)) set.seed(seed)
    epsgen <- replicate(R, do.call(errorgen, c("n" = n, errorparam)),
                        simplify = FALSE)
    egen <- lapply(epsgen, function(eps) M %*% eps)

    statgen <- vapply(1:R, function(j) do.call(what = hettest,
                args = append(list("mainlm" = list("X" = X, "e" = egen[[j]]),
                         "statonly" = TRUE), arguments)), NA_real_)
  } else { # Goldfeld-Quandt
    if (exists("method", where = environment(), inherits = FALSE)) {
      method <- match.arg(method, c("parametric", "nonparametric"))
      if (method == "nonparametric") stop("This method is only available for tests with a continuous test statistic")
    }

    if (!is.null(seed)) set.seed(seed)
    epsgen <- replicate(R, do.call(errorgen, c("n" = n, errorparam)),
                        simplify = FALSE)

    hasintercept <- columnof1s(X)
    if (class(mainlm) == "list") {
      if (hasintercept[[1]]) {
        if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
        colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
      } else {
        colnames(X) <- paste0("X", 1:p)
      }
    }

    if (!exists("deflator", where = environment(), inherits = FALSE)) deflator <- NULL
    if (!exists("prop_central", where = environment(), inherits = FALSE)) prop_central <- 1 / 3
    if (!exists("group1prop", where = environment(), inherits = FALSE)) group1prop <- 1 / 2
    checkdeflator(deflator, X, p, hasintercept[[1]])
    theind <- gqind(n, prop_central, group1prop)

    if (!is.null(deflator)) {
      X <- X[order(X[, deflator]), , drop = FALSE]
    }
    M1 <- fastM(X[theind[[1]], , drop = FALSE], length(theind[[1]]))
    M2 <- fastM(X[theind[[2]], , drop = FALSE], length(theind[[2]]))
    thedf2 <- (length(theind[[2]]) - p)
    thedf1 <- (length(theind[[1]]) - p)

    egen1 <- lapply(epsgen, function(eps) M1 %*% eps[theind[[1]]])
    egen2 <- lapply(epsgen, function(eps) M2 %*% eps[theind[[2]]])

    statgen <- unlist(mapply(FUN = function(e1, e2)
      (sum(e2 ^ 2) / thedf2) / (sum(e1 ^ 2) / thedf1), egen1, egen2,
      SIMPLIFY = FALSE))
  }

  if (alternative == "greater") {
    teststat <- sum(statgen >= statobs)
  } else if (alternative == "less") {
    teststat <- sum(statgen <= statobs)
  } else if (alternative == "two.sided") {
    teststat <- min(sum(statgen >= statobs),
                    sum(statgen <= statobs))
  }

  pval <- (teststat + 1) / (R + 1) *
    ifelse(alternative == "two.sided", 2, 1)

  rval <- structure(list(statistic = teststat, parameter = R, p.value = pval,
               null.value = "Homoskedasticity",
               alternative = alternative, method = "Monte Carlo"),
               class = "htest")
  broom::tidy(rval)
}
