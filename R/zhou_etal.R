#' Zhou, Song, and Thompson's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the methods of
#'    \insertCite{Zhou15;textual}{skedastic} for testing for
#'    heteroskedasticity in a linear regression model. The authors propose
#'    three variations based on comparisons between sandwich and model-based
#'    estimators for the variances of individual regression coefficient
#'    esimators. The covariate-specific method computes a test statistic and
#'    \eqn{p}-value for each column of the auxiliary design matrix (which is,
#'    by default, the original design matrix with intercept omitted). The
#'    \eqn{p}-values undergo a Bonferroni correction to control overall test
#'    size. When the null hypothesis is rejected in this case, it also provides
#'    information about which auxiliary design variable is associated with the
#'    error variance. The pooled method computes a single test statistic and
#'    \eqn{p}-value and is thus an omnibus test. The hybrid method returns the
#'    minimum \eqn{p}-value between the Bonferroni-corrected covariate-specific
#'    \eqn{p}-values and the pooled \eqn{p}-value, multiplying it by 2 for a
#'    further Bonferroni correction. The test statistic returned is that
#'    which corresponds to the minimum \eqn{p}-value. The covariate-specific
#'    and pooled test statistics have null distributions that are
#'    asymptotically normal with mean 0. However, the variance is intractable
#'    and thus perturbation sampling is used to compute \eqn{p}-values
#'    empirically.
#'
#' @param method A character specifying which of the three test methods to
#'    implement; one of \code{"pooled"}, \code{"covariate-specific"}, or
#'    \code{"hybrid"} (which combines the other two). Partial matching is
#'    used.
#' @param Bperturbed An integer specifying the number of perturbation samples
#'    to generate when estimating the \eqn{p}-value. Defaults to \code{500L}.
#' @param seed An integer specifying a seed to pass to
#'    \code{\link[base]{set.seed}} for random number generation. This allows
#'    for reproducibility of perturbation sampling. A value of \code{NA}
#'    results in not setting a seed.
#'
#' @inheritParams breusch_pagan
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' zhou_etal(mtcars_lm, method = "pooled")
#' zhou_etal(mtcars_lm, method = "covariate-specific")
#' zhou_etal(mtcars_lm, method = "hybrid")
#'

zhou_etal <- function(mainlm, auxdesign = NA,
                      method = c("pooled", "covariate-specific", "hybrid"),
                      Bperturbed = 500L, seed = 1234, statonly = FALSE) {

  method <- match.arg(method, c("pooled", "covariate-specific", "hybrid"))

  auxfitvals <- ifelse(all(is.na(auxdesign)) | is.null(auxdesign), FALSE,
                                    auxdesign == "fitted.values")
  processmainlm(m = mainlm, needy = auxfitvals, needyhat = auxfitvals,
                needp = FALSE)

  if (all(is.na(auxdesign)) || is.null(auxdesign)) {
    Z <- X
  } else if (is.character(auxdesign)) {
    if (auxdesign == "fitted.values") {
      Z <- t(t(yhat))
    } else stop("Invalid character value for `auxdesign`")
  } else {
    Z <- auxdesign
    if (nrow(auxdesign) != nrow(X)) stop("No. of observations in `auxdesign`
                                         must match\nno. of observations in
                                         original model.")
  }

  hasintercept <- columnof1s(Z)
  if (hasintercept[[1]]) {
    Z <- Z[, -hasintercept[[2]], drop = FALSE]
  }

  q <- ncol(Z)
  if (q == 1 && method %in% c("covariate-specific", "hybrid")) {
    message("Auxiliary design matrix consists of only one column; method changed to `pooled`")
    method <- "pooled"
  }
  n <- nrow(X)

  # y <- y - mean(y)
  # X <- apply(X, 2, function(x) x - mean(x))
  # newlm <- stats::lm.fit(X, y)
  # e <- newlm$residuals

  if (!is.na(seed)) set.seed(seed)
  Xi <- replicate(Bperturbed, stats::rnorm(n), simplify = FALSE)
  sigmahatsq <- sum(e ^ 2) / n
  H <- Z %*% solve(crossprod(Z)) %*% t(Z)

  if (method == "pooled" || method == "hybrid") {
    wpool <- diag(H) / q
    IRpool <- sum(e ^ 2 * wpool) / sigmahatsq
    Wpoolobs <- sqrt(n) * (IRpool - 1)
    if (statonly && method == "pooled") return(Wpoolobs)
    Wpoolstar <- vapply(Xi, function(xi) sqrt(n) * sum(((wpool - IRpool / n) *
                  (e ^ 2 / sigmahatsq - 1) - (IRpool - 1) / n) * xi), NA_real_)
    pvalpool <- sum(Wpoolstar >= Wpoolobs) / Bperturbed
  }
  if (method == "covariate-specific" || method == "hybrid") {

    Hminus <- lapply(1:q, function(j) Z[, -j, drop = FALSE] %*%
              solve(crossprod(Z[, -j, drop = FALSE])) %*%
                t(Z[, -j, drop = FALSE]))
    w <- lapply(1:q, function(j) diag(H) - diag(Hminus[[j]]))
    IR <- vapply(1:q, function(j) sum(e ^ 2 * w[[j]]) / sigmahatsq, NA_real_)
    Wobs <- sqrt(n) * (IR - 1)
    if (statonly && method == "covariate-specific") return(Wobs)
    Wstar <- lapply(1:q, function(j) vapply(Xi,
              function(xi) sqrt(n) * sum(((w[[j]] - IR[j] / n) *
               (e ^ 2 / sigmahatsq - 1) - (IR[j] - 1) / n) * xi), NA_real_))
    pvalcov <- q * vapply(1:q, function(j) sum(Wstar[[j]] >= Wobs[j]) /
                 Bperturbed, NA_real_)
    pvalcov[pvalcov > 1] <- 1
    message("covariate-specific p-values include Bonferroni correction")
  }

  if (method == "pooled") {
    teststat <- Wpoolobs
    pval <- pvalpool
  } else if (method == "covariate-specific") {
    teststat <- Wobs
    pval <- pvalcov
  } else if (method == "hybrid") {
    minpval <- min(pvalpool, pvalcov)
    pval <- ifelse(2 * minpval > 1, 1, 2 * minpval)
    message("hybrid p-value includes Bonferroni correction")
    if (minpval == pvalpool) {
      teststat <- Wpoolobs
    } else {
      teststat <- Wobs
    }
    if (statonly) return(teststat)
  }

  rval <- structure(list(statistic = teststat, parameter = NULL,
               p.value = pval,
               null.value = "Homoskedasticity",
               alternative = "Heteroskedasticity", method = method), class = "htest")
  broom::tidy(rval)
}
