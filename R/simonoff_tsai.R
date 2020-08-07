#' Simonoff-Tsai Tests for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the modified profile likelihood ratio test and
#'   score test of \insertCite{Simonoff94;textual}{skedastic} for testing
#'   for heteroskedasticity in a linear regression model.
#'
#' The Simonoff-Tsai Likelihood Ratio Test involves a modification of the
#' profile likelihood function so that the nuisance parameter will be
#' orthogonal to the parameter of interest. The maximum likelihood estimate
#' of \eqn{\lambda} (called \eqn{\delta} in
#' \insertCite{Simonoff94;textual}{skedastic}) is computed from the modified
#' profile log-likelihood function using the Nelder-Mead algorithm in
#' \code{\link[maxLik]{maxLik}}. Under the null hypothesis of homoskedasticity,
#' the distribution of the test statistic is asymptotically chi-squared with
#' \eqn{q} degrees of freedom. The test is right-tailed.
#'
#' The Simonoff-Tsai Score Test entails adding a term to either the score
#' statistic of \insertCite{Cook83;textual}{skedastic} (a test implemented in
#' \code{\link{cook_weisberg}}) or to that of
#' \insertCite{Koenker81;textual}{skedastic} (a test implemented in
#' \code{\link{breusch_pagan}} with argument \code{koenker} set to
#' \code{TRUE}), in order to improve the robustness of these respective tests
#' in the presence of non-normality. This test likewise has a test statistic
#' that is asymptotically \eqn{\chi^2(q)}-distributed and the test is likewise
#' right-tailed.
#'
#' The assumption of underlying both tests is that
#' \eqn{\mathrm{Cov}(\epsilon)=\sigma^2 W}, where \eqn{W} is
#' an \eqn{n\times n} diagonal matrix with \eqn{i}th diagonal element
#' \eqn{w_i=w(Z_i, \lambda)}. Here, \eqn{Z_i} is the \eqn{i}th row of an
#' \eqn{n \times q} nonstochastic auxiliary design matrix \eqn{Z}. Note:
#' \eqn{Z} as defined here does not have a column of ones, but is concatenated
#' to a column of ones when used in an auxiliary regression.
#' \eqn{\lambda} is a \eqn{q}-vector of unknown parameters, and \eqn{w(\cdot)}
#' is a real-valued, twice-differentiable function having the property that
#' there exists some \eqn{\lambda_0} for which
#' \eqn{w(Z_i,\lambda_0)=0} for all \eqn{i=1,2,\ldots,n}. Thus, the null
#' hypothesis of homoskedasticity may be expressed as
#' \eqn{\lambda=\lambda_0}.
#'
#' In the score test, the added term in the test statistic is of the
#' form
#' \deqn{\sum_{j=1}^{q} \left(\sum_{i=1}^{n} h_{ii} t_{ij}\right) \tau_j},
#' where \eqn{t_{ij}} is the \eqn{(i,j)}th element of the Jacobian matrix
#' \eqn{J} evaluated at \eqn{\lambda=\lambda_0}:
#' \deqn{t_{ij}=\left.\frac{\partial w(Z_i, \lambda)}{\partial \lambda_j}\right|_{\lambda=\lambda_0}},
#' and \eqn{\tau=(\bar{J}'\bar{J})^{-1}\bar{J}'d}, where \eqn{d} is the
#' \eqn{n}-vector whose \eqn{i}th element is \eqn{e_i^2\bar{\sigma}^{-2}},
#' \eqn{\bar{\sigma}^2=n^{-1}e'e}, and \eqn{\bar{J}=(I_n-1_n 1_n'/n)J}.
#'
#' @param method A character specifying which of the tests proposed in
#'    \insertCite{Simonoff94;textual}{skedastic} to implement. \code{"mlr"}
#'    corresponds to the modified profile likelihood ratio test, and
#'    \code{"score"} corresponds to the score test.
#' @param basetest A character specifying the base test statistic which is
#'    robustified using the added term described in Details. \code{"koenker"}
#'    corresponds to the test statistic produced by \code{\link{breusch_pagan}}
#'    with argument \code{koenker} set to \code{TRUE}, while
#'    \code{"cook_weisberg"} corresponds to the test statistic produced by
#'    \code{\link{cook_weisberg}}. Partial matching is used. This argument is
#'    only used if \code{method} is \code{"score"}.
#' @param bartlett A logical specifying whether a Bartlett correction should be
#'    made, as per \insertCite{Ferrari04;textual}{skedastic}, to improve the
#'    fit of the test statistic to the asymptotic null distribution. This
#'    argument is only applicable where \code{method} is \code{"mlr"}, and is
#'    implemented only where \code{hetfun} is \code{"mult"} or
#'    \code{"logmult"}.
#' @param optmethod A character specifying the optimisation method to use with
#'    \code{\link[stats]{optim}}, if \code{method} is \code{"mlr"}. The
#'    default, \code{"Nelder-Mead"}, corresponds to the default \code{method}
#'    value in \code{\link[stats]{optim}}. Warnings about Nelder-Mead algorithm
#'    being unreliable for one-dimensional optimization have been suppressed,
#'    since the algorithm does appear to work for the three implemented choices
#'    of \code{hetfun}.
#' @param ... Optional arguments to pass to \code{\link[stats]{optim}}, such as
#'    \code{par} (initial value of \eqn{\lambda}) and \code{maxit} (maximum
#'    number of iterations to use in optimisation algorithm), and \code{trace}
#'    (to provide detailed output on optimisation algorithm). Default initial
#'    value of \eqn{\lambda} is \code{rep(1e-3, q)}.
#' @inheritParams breusch_pagan
#' @inheritParams cook_weisberg
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
#' simonoff_tsai(mtcars_lm, method = "score")
#' simonoff_tsai(mtcars_lm, method = "score", basetest = "cook_weisberg")
#' simonoff_tsai(mtcars_lm, method = "mlr")
#' simonoff_tsai(mtcars_lm, method = "mlr", bartlett = FALSE)
#' \dontrun{simonoff_tsai(mtcars_lm, auxdesign = data.frame(mtcars$wt, mtcars$qsec),
#'  method = "mlr", hetfun = "logmult")}
#'

simonoff_tsai <- function(mainlm, auxdesign = NA, method = c("mlr", "score"),
                  hetfun = c("mult", "add", "logmult"),
                  basetest = c("koenker", "cook_weisberg"), bartlett = TRUE,
                  optmethod = "Nelder-Mead", statonly = FALSE,
                  ...) {

  basetest <- match.arg(basetest, c("koenker", "cook_weisberg"))
  method <- match.arg(method, c("mlr", "score"))
  hetfun <- match.arg(hetfun, c("mult", "add", "logmult"))

  auxfitvals <- ifelse(all(is.na(auxdesign)) | is.null(auxdesign), FALSE,
                                    auxdesign == "fitted.values")
  processmainlm(m = mainlm, needy = (method == "mlr"), needyhat = auxfitvals,
                needp = TRUE)

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
  n <- nrow(Z)

  if (method == "score") {
    sigma_barsq <- sum(e ^ 2) / n
    d <- e ^ 2 / sigma_barsq

    if (hetfun == "mult") {
      J <- Z
    } else if (hetfun == "logmult") {
      if (any(Z <= 0)) stop("When `hetfun` is `\"logmult\"`, auxiliary design matrix elements must be strictly positive")
      J <- log(Z)
    } else if (hetfun == "add") {
      J <- 2 * Z
    } else stop("Invalid hetfun argument")

    if (basetest == "koenker") {

      w_hat <- e ^ 2 - sigma_barsq
      basestat <- n * sum(stats::lm.fit(cbind(1, Z), w_hat)$fitted.values ^ 2) /
        sum(w_hat ^ 2)
    } else if (basetest == "cook_weisberg") {
      auxres <- stats::lm.fit(cbind(1, J), d)$residuals
      basestat <- (sum(d ^ 2) - n * mean(d) ^ 2
                   - sum(auxres ^ 2)) / 2
    }

    tij <- sapply(1:q, function(j) vapply(1:n, function(i) J[i, j] -
                                            sum(J[, j]) / n, NA_real_), simplify = TRUE)
    Jbar <- (diag(n) - matrix(data = 1 / n, nrow = n, ncol = n)) %*% J
    tau <- solve(crossprod(Jbar)) %*% t(Jbar) %*% d
    h <- diag(X %*% solve(crossprod(X)) %*% t(X))
    addterm <- sum(vapply(1:q, function(j) sum(h * tij[, j]) * tau[j], NA_real_))
    teststat <- basestat + addterm
    if (statonly) return(teststat)
  } else if (method == "mlr") {

    if (hetfun == "mult") {
      w <- function(Zi, lambda) {
        force(exp(sum(lambda * Zi)))
      }
    } else if (hetfun == "logmult") {
      w <- function(Zi, lambda) {
        if (any(Zi <= 0)) stop("Zi cannot be negative")
        force(exp(sum(lambda * log(Zi))))
      }
    } else if (hetfun == "add") {
      w <- function(Zi, lambda) {
        force((1 + sum(lambda * Zi)) ^ 2)
      }
    } else stop("Invalid hetfun argument")
    ellp <- function(lambda) {
      W <- diag(vapply(1:n, function(i) w(Z[i, ], lambda), NA_real_))
      Winv <- solve(W)
      beta <- solve(t(X) %*% Winv %*% X) %*% t(X) %*% Winv %*% y
      sigmasq <- t(y - X %*% beta) %*% Winv %*% (y - X %*% beta) / n
      (- n / 2 * log(sigmasq) - 1 / 2 * sum(vapply(1:n, function(i)
        log(w(Z[i, ], lambda)), NA_real_)) - (2 * sigmasq) ^ (-1) *
        t(y - X %*% beta) %*% Winv %*% (y - X %*% beta))
    }
    lambda0 <- rep(0, q)
    lambda_start <- rep(1e-3, q)

    arguments <- list(...)
    invisible(list2env(arguments, envir = environment()))

    if (exists("par", where = environment(), inherits = FALSE)) {
      if (length(par) == 1) {
        lambda_start <- rep(par, q)
      } else {
        lambda_start <- par
      }
      arguments$par <- NULL
    }

    controllist <- append(list("fnscale" = -1, "warn.1d.NelderMead" = FALSE),
                          arguments)

    MLE <- stats::optim(par = lambda_start, fn = ellp, method = optmethod,
                        control = controllist, hessian = TRUE)

    if (MLE$convergence != 0) warning("Algorithm for likelihood maximisation did not converge; consider trying a larger maxit value")

    # constraints = list("ineqA" = diag(q), "ineqB" = rep(0, q))
    L <- -2 * (ellp(lambda0) - MLE$value)
    w_at_MLE <- vapply(1:n, function(i) w(Z[i, ], MLE$par), NA_real_)
    Ghat <- diag(w_at_MLE / (prod(w_at_MLE ^ (1 / n))))
    Xmhat <- pracma::sqrtm(solve(Ghat)) %*% X
    teststat <- as.double((n - p - 2) / n * L +
      log(det(crossprod(X)) / det(crossprod(Xmhat))))
    if (bartlett) {
      if (hetfun == "add") {
        warning("Bartlett correction not implemented for this choice of hetfun")
        thecm <- 0
      } else {
        if (hetfun == "mult") {
          Zeta <- as.matrix(Z)
        } else if (hetfun == "logmult") {
          if (any(Z <= 0)) stop("All elements of auxiliary design matrix must be positive for this choice of hetfun")
          Zeta <- as.matrix(log(Z))
        }
        Zetabar <- colMeans(Zeta)
        Hferr <- (Zeta - Zetabar) %*% solve(crossprod(Zeta -
                         Zetabar)) %*% t(Zeta - Zetabar)
        thecm <- -1/2 * sum(diag(Hferr) ^ 2) + 1/(2*n) * (sum(diag(Hferr)) ^ 2) +
          1/2 * sum(vapply(1:n, function(l) Hferr[l, l] *
            sum(vapply(1:n, function(m) Hferr[l, m] * Hferr[m, m], NA_real_)),
            NA_real_)) + 1/3 * sum(Hferr ^ 3) - 2 / n * sum(diag(Hferr)) +
            1 / n * sum(Hferr ^ 2)
      }
      teststat <- teststat / (1 + thecm / q)
    }
    if (statonly) return(teststat)
  }

  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = method),
                         class = "htest")
  broom::tidy(rval)
}
