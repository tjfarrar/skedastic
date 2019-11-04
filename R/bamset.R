#' Ramsey's BAMSET Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the Bartlett's \eqn{M} Specification Error Test
#'    (BAMSET) method of \insertCite{Ramsey69;textual}{skedastic} for testing
#'    for heteroskedasticity in a linear regression model.
#'
#' BAMSET is an analogue of Bartlett's \eqn{M} Test for heterogeneity of
#'    variances across independent samples from \eqn{k} populations. In this
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
#'    argument must be specified and cannot be left as \code{NULL}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams goldfeld_quandt
#' @inheritParams blus
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' n <- 20
#' p <- 4
#' set.seed(9586)
#' X <- matrix(data = runif(n * (p - 1)), nrow = n, ncol = p - 1)
#' # Response values generated under homoskedasticity
#' y_H0 <- rnorm(n, mean = 1 + rowSums(X), sd = 1)
#' bamset(lm(y_H0 ~ X))
#'# Response values generated under heteroskedasticity associated with X
#' y_HA <- rnorm(n, mean = 1 + rowSums(X), sd = X[, 1] ^ 2)
#' bamset(lm(y_HA ~ X))
#' bamset(lm(y_HA ~ X), deflator = "X1")
#'

bamset <- function(mainlm, k = 3, deflator = NULL, correct = TRUE,
                   omitatmargins = TRUE, omit = NULL) {

  if (!is.null(deflator)) {
    if (class(mainlm) == "lm") {
      y <- model.response(model.frame(mainlm))
      X <- model.matrix(mainlm)
      p <- ncol(X)
      n <- length(y)
      hasintercept <- columnof1s(X)
    } else if (class(mainlm) == "list") {
      y <- mainlm[[1]]
      X <- mainlm[[2]]
      badrows <- which(apply(cbind(mainlm[[1]], mainlm[[2]]), 1, function(x)
        any(is.na(x), is.nan(x), is.infinite(x))))
      if (length(badrows) > 0) {
        stop("NA/NaN/Inf values not permitted in data")
      }
      p <- ncol(X)
      n <- length(y)
      hasintercept <- columnof1s(X)
      if (hasintercept[[1]]) {
        if (hasintercept[[2]] != 1) stop("Column of 1's must be first column
                                         of design matrix")
        colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
      } else {
        colnames(X) <- paste0("X", 1:p)
      }
    }

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
    } else stop("`deflator` must be integer or character")

    y <- y[order(X[, deflator])]
    X <- X[order(X[, deflator]), ]
    mainlm <- list(y, X)

  } else {
    if (class(mainlm) == "lm") {
      n <- length(mainlm$residuals)
      p <- length(mainlm$coefficients)
    } else if (class(mainlm) == "list") {
      n <- length(mainlm[[1]])
      p <- ncol(mainlm[[2]])
      badrows <- which(apply(cbind(mainlm[[1]], mainlm[[2]]), 1, function(x)
        any(is.na(x), is.nan(x), is.infinite(x))))
      if (length(badrows) > 0) {
        stop("NA/NaN/Inf values not permitted in data")
      }
    }
  }

  nprime <- n - p
  min_subset_size <- as.integer(nprime / k) # called r_1 in Ramsey (1969)
  nprime_modk <- nprime %% k

  # slightly different definition of v_i compared with Ramsey (1969)
  # this one makes the subsets more equitable in size (not differing by more than 1)
  v <- rep(min_subset_size, k) + c(rep(1, nprime_modk), rep(0, k - nprime_modk))

  if (omitatmargins) omit <- margin_indices(v, p)
  res <- blus(mainlm, omit)

  s_sq <- unlist(lapply(sub_ind, function(i, e)
    sum(e[i] ^ 2, na.rm = T), e = res)) / v
  s_sq_tot <- sum(res ^ 2, na.rm = T) / nprime

  teststat <- nprime * log(s_sq_tot) - sum(v * log(s_sq))
  if (correct) {
    scaling_constant <- 1 + sum(1 / v - 1 / nprime) / (3 * (k - 1))
    teststat <- teststat / scaling_constant
  }

  df <- k - 1
  pval <- 1 - pchisq(teststat, df = df)

  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "Heteroskedasticity", method = "BAMSET"),
                    class = "htest")
  broom::tidy(rval)
}
