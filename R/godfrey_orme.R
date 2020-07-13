#' Godfrey and Orme's Nonparametric Bootstrap Test for Heteroskedasticity in a Linear
#' Regression Model
#'
#' This function implements the method of
#'    \insertCite{Godfrey99;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model. The procedure is more
#'    clearly described in \insertCite{Godfrey06;textual}{skedastic}.
#'
#' The procedure runs as follows. (1) The observed
#'    value of the test statistic \eqn{T_0} is computed using \code{hettest}.
#'    (2) A sample \eqn{e_1^\star,e_2^\star,\ldots,e_n^\star} is drawn with
#'    replacement from the OLS residuals. (3) Bootstrapped response values are
#'    computed as \eqn{y_i^{\star}=x_i' \hat{\beta}+e_i^\star,i=1,2,\ldots,n}.
#'    (4) Bootstrapped test statistic value \eqn{T^\star} is computed from the
#'    regression of \eqn{y^\star} on \eqn{X} using function \code{hettest}.
#'    (5) Steps (2)-(4) are repeated until \eqn{B} bootstrapped test statistic
#'    values are computed. (6) Empirical \eqn{p}-value is computed by comparing
#'    the bootstrapped test statistic values to the observed test statistic
#'    value. Note that, if \code{alternative} is set to \code{"two.sided"}, the
#'    one-sided \eqn{p}-value is doubled (\code{\link{twosidedpval}} cannot
#'    be used in this case).
#'
#' @param B An integer specifying the number of nonparametric bootstrap samples
#'    to generate. Defaults to \code{1000L}.
#' @param hettest A function, or a character specifying the name of a function,
#'    that implements a heteroskedasticity test on a linear regression model.
#'    The function is called with the \code{statonly} argument set to
#'    \code{TRUE} to improve computational efficiency.
#' @param seed An integer specifying a seed to pass to
#'    \code{\link[base]{set.seed}} for random number generation. This allows
#'    reproducibility of bootstrap results. A value of \code{NULL}
#'    results in not setting a seed.
#' @param ... Additional arguments to pass to \code{hettest}
#'
#' @inheritParams breusch_pagan
#' @inheritParams dufour_etal
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object
#'    is not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' godfrey_orme(mtcars_lm, hettest = breusch_pagan)
#'

godfrey_orme <- function(mainlm, hettest, B = 1000L, alternative = c("greater",
                        "less", "two.sided"), seed = 1234, ...) {

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  if (is.character(hettest)) hettest <- get(hettest)
  arguments <- list(...)
  invisible(list2env(arguments, envir = environment()))
  if ("alternative" %in% names(formals(hettest))) arguments$alternative <- alternative
  hettestchar <- deparse(substitute(hettest))
  processmainlm(m = mainlm)
  n <- length(e)

  if (exists("deflator", where = environment(), inherits = FALSE)) {
    if (is.character(deflator))
      arguments$deflator <- which(colnames(X) == deflator)
  }

  statobs <- do.call(what = hettest,
              args = append(list("mainlm" = list("y" = y, "X" = X),
                                 "statonly" = TRUE), arguments))

  if (!is.null(seed)) set.seed(seed)
  estar <- replicate(B, sample(e, size = n, replace = TRUE), simplify = FALSE)
  betahat <- solve(crossprod(X)) %*% t(X) %*% y
  ystar <- lapply(estar, function(u) X %*% betahat + u)

  statgen <- vapply(1:B, function(b) do.call(what = hettest,
              args = append(list("mainlm" = list("X" = X, "y" = ystar[[b]]),
                       "statonly" = TRUE), arguments)), NA_real_)

  if (alternative == "greater") {
    teststat <- sum(statgen >= statobs)
  } else if (alternative == "less") {
    teststat <- sum(statgen <= statobs)
  } else if (alternative == "two.sided") {
    teststat <- min(sum(statgen >= statobs),
                    sum(statgen <= statobs))
  }

  pval <- teststat / B * ifelse(alternative == "two.sided", 2, 1)

  rval <- structure(list(statistic = teststat, parameter = B, p.value = pval,
               null.value = "Homoskedasticity",
               alternative = alternative, method = "Nonpar. Bootstrap"),
               class = "htest")
  broom::tidy(rval)
}
