#' Ramsey's BAMSET Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the Bartlett's \eqn{M} Specification Error Test
#'    (BAMSET) method of \insertCite{Ramsey69;textual}{skedastic} for testing
#'    for heteroskedasticity in a linear regression model.
#'
#' @details BAMSET is an analogue of Bartlett's \eqn{M} Test for heterogeneity
#'    of variances across independent samples from \eqn{k} populations. In this
#'    case the populations are \eqn{k} subsets of the residuals from a linear
#'    regression model. In order to meet the independence assumption,
#'    \link[=blus]{BLUS residuals} are computed, meaning that only \eqn{n-p}
#'    observations are used (where \eqn{n} is the number of rows and \eqn{p}
#'    the number of columns in the design matrix). Under the null hypothesis
#'    of homoskedasticity, the test statistic is asymptotically chi-squared
#'    distributed with \eqn{k-1} degrees of freedom. The test is right-tailed.
#'
#' @param k An integer. The number of subsets (>= 2) into which the BLUS residuals are
#'    to be partitioned. Defaults to 3, the value suggested in
#'    \insertCite{Ramsey69;textual}{skedastic}.
#' @param correct A logical. Should the test statistic be divided by a scaling
#'    constant to improve the chi-squared approximation? Defaults to
#'    \code{TRUE}.
#' @param omitatmargins A logical. Should the indices of observations at the
#'    margins of the \code{k} subsets be passed to \code{\link{blus}} as the
#'    \code{omit} argument? If \code{TRUE} (the default), this overrides any
#'    \code{omit} argument passed directly. If \code{FALSE}, the \code{omit}
#'    argument must be specified and cannot be left as \code{NA}. If
#'    \code{categorical} is \code{TRUE}, setting \code{omitatmargins} to
#'    \code{TRUE} results in omitting observations from the most frequently
#'    occurring factor levels or values.
#' @param categorical A logical. Is the deflator a categorical variable? If
#'    so, the number of levels will be used as \eqn{k} with each level forming
#'    a subset. Defaults to \code{FALSE}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams goldfeld_quandt
#' @inheritParams blus
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
#' bamset(mtcars_lm, deflator = "qsec", k = 3)
#'
#' # BLUS residuals cannot be computed with given `omit` argument and so
#' # omitted indices are randomised:
#' bamset(mtcars_lm, deflator = "qsec", k = 4, omitatmargins = FALSE, omit = "last")
#'

bamset <- function(mainlm, k = 3, deflator = NA, correct = TRUE,
                   omitatmargins = TRUE, omit = NA,
                   categorical = FALSE, statonly = FALSE) {

  processmainlm(m = mainlm, needy = FALSE)

  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column
                                         of design matrix")
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

  n <- nrow(X)
  nprime <- n - p

  if (categorical) {
    if (is.factor(X[, deflator])) {
      lev <- levels(X[, deflator])
    } else {
      lev <- unique(sort(X[, deflator]))
    }
    k <- length(lev)
    min_subset_size <- as.integer(nprime / k) # called r_1 in Ramsey (1969)
    nprime_modk <- nprime %% k
    v <- table(X[, deflator])

  } else {
    min_subset_size <- as.integer(nprime / k) # called r_1 in Ramsey (1969)
    nprime_modk <- nprime %% k

    # slightly different definition of v_i compared with Ramsey (1969)
    # this one makes the subsets more equitable in size (not differing by more than 1)
    v <- rep(min_subset_size, k) + c(rep(1, nprime_modk), rep(0, k - nprime_modk))
  }
  sub_ind <- vector("list", k)
  sub_ind[[1]] <- seq.int(from = 1, by = 1, length.out = v[1])
  for (j in 2:k) {
    sub_ind[[j]] <- seq.int(from = max(sub_ind[[j - 1]]) + 1, by = 1,
                            length.out = v[j])
  }

  if (omitatmargins)
    omit <- margin_indices(v, p, sub_ind,
                                       categ = categorical)

  res <- blus(mainlm = list("X" = X, "e" = e), omit)

  s_sq <- unlist(lapply(sub_ind, function(i, e)
    sum(e[i] ^ 2, na.rm = TRUE), e = res)) / v
  s_sq_tot <- sum(res ^ 2, na.rm = TRUE) / nprime

  teststat <- nprime * log(s_sq_tot) - sum(v * log(s_sq))
  if (correct) {
    scaling_constant <- 1 + (sum(1 / v) - 1 / nprime) / (3 * (k - 1))
    teststat <- teststat / scaling_constant
  }

  if (statonly) return(teststat)

  df <- k - 1
  pval <- stats::pchisq(teststat, df = df, lower.tail = FALSE)

  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = "BAMSET"),
                    class = "htest")
  broom::tidy(rval)
}
