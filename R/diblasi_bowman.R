#' Diblasi and Bowman's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the nonparametric test of
#'    \insertCite{Diblasi97;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' The test entails undertaking a transformation of the OLS residuals
#'   \eqn{s_i=\sqrt{|e_i|}-E_0(\sqrt{|e_i|})}, where \eqn{E_0} denotes
#'   expectation under the null hypothesis of homoskedasticity. The kernel
#'   method of nonparametric regression is used to fit the relationship
#'   between these \eqn{s_i} and the explanatory variables. This leads to a
#'   test statistic \eqn{T} that is a ratio of quadratic forms involving the
#'   vector of \eqn{s_i} and the matrix of normal kernel weights. Although
#'   nonparametric in its method of fitting the possible heteroskedastic
#'   relationship, the distributional approximation used to compute
#'   \eqn{p}-values assumes normality of the errors.
#'
#' @param H A hyperparameter denoting the bandwidth matrix in the kernel
#'   function used for weights in nonparametric smoothing. If a double of
#'   length 1 (the default), \code{H} is set to \eqn{h I_{p^\prime}} where
#'   \eqn{h} is the scalar bandwidth value entered and \eqn{I_{p^\prime}}
#'   is the \eqn{p^prime \times p^\prime} identity matrix (where
#'   \eqn{p^\prime} is the number of columns in the \eqn{X} matrix, excluding
#'   an intercept if present). If a double of length \eqn{p^\prime}, \code{H}
#'   is set to \eqn{diag(h)} where \eqn{h} is the bandwidth vector entered.
#'   If \code{H} is a \eqn{p^\prime\times p^\prime} matrix it is used as is.
#'   Any other dimensionality of \code{H} results in an error.
#' @param distmethod A character specifying the method by which to estimate
#'   the \eqn{p}-values, either \code{"moment.match"} or \code{"bootstrap"}.
#' @param B An integer specifying the number of nonparametric bootstrap
#'   replications to be used, if \code{distmethod="bootstrap"}.
#' @param ignorecov A logical. If \code{TRUE} (the default), the
#'   variance-covariance matrix of \eqn{s} is assumed to be diagonal. (This
#'   assumption is, strictly speaking, invalid, but improves computational
#'   efficiency considerably for what should be a small loss of
#'   statistical precision).
#'
#' @inheritParams breusch_pagan
#' @inheritParams wilcox_keselman
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
#' diblasi_bowman(mtcars_lm)
#' diblasi_bowman(mtcars_lm, ignorecov = FALSE)
#' diblasi_bowman(mtcars_lm, distmethod = "bootstrap")
#'
#' # Example discussed in Diblasi and Bowman (1997)
#' malecats_lm <- lm(Hwt ~ Bwt, data = boot::catsM)
#' diblasi_bowman(malecats_lm, H = (max(boot::catsM$Bwt) - min(boot::catsM$Bwt)) / 8)
#'

diblasi_bowman <- function(mainlm, distmethod = c("moment.match", "bootstrap"),
                            H = 0.08, ignorecov = TRUE, B = 500L, seed = 1234,
                           statonly = FALSE) {

  distmethod <- match.arg(distmethod, c("moment.match", "bootstrap"))
  processmainlm(m = mainlm, needy = (distmethod == "bootstrap"),
                needyhat = (distmethod == "bootstrap"))

  hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) {
    if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
    colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
  } else {
    colnames(X) <- paste0("X", 1:p)
  }
  n <- nrow(X)
  H <- as.matrix(H)
  if (max(dim(H)) == 1) {
    H <- as.double(H)
  } else if (min(dim(H)) == 1) {
    if (length(as.double(H)) != (p - hasintercept[[1]])) {
      stop("Invalid dimensionality of H")
    }
    H <- diag(as.double(H))
  } else if (min(dim(H)) > 1) {
    if (any(dim(H) != (p - hasintercept[[1]]))) {
      stop("Invalid dimensionality of H")
    }
  }

  if (hasintercept[[1]]) {
    Xforweights <- X[, -1, drop = FALSE]
  } else {
    Xforweights <- as.matrix(X)
  }

  if (is.matrix(H)) {
    Hinv <- solve(H)
    w <- function(xi, xj) exp(-1 / 2 * (t(xi - xj) %*% Hinv %*% Hinv %*% (xi - xj)))
  } else {
    w <- function(xi, xj) exp(-1 / (2 * H ^ 2) * (t(xi - xj) %*% (xi - xj)))
  }
  W <- apply(Xforweights, 1, function(xi) apply(Xforweights, 1, function(xj) w(xi, xj)))
  W <- W / rowSums(W)
  Bmat <- (t(diag(n) - W)) %*% (diag(n) - W)
  Cmat <- diag(n) - matrix(data = 1 / n, nrow = n, ncol = n) - Bmat
  sigma_hat_sq <- sum(e ^ 2) / (n - p)
  M <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
  E0 <- gamma(3 / 4) * (2 * sigma_hat_sq * diag(M)) ^ (1 / 4) / sqrt(pi)
  s <- sqrt(abs(e)) - E0
  teststat <- as.double((t(s) %*% Cmat %*% s) / (t(s) %*% Bmat %*% s))
  if (statonly) return(teststat)

  if (distmethod == "moment.match") {
    if (ignorecov) {
      Sigma <- diag(sqrt(2) / pi * (sqrt(pi) - gamma(3 / 4) ^ 2) *
        sqrt(sigma_hat_sq * diag(M)))
    } else {
      # Sigma <- matrix(data = unlist(lapply(1:n, function(i) lapply(1:n,
      #   function(j) ifelse(j <= i, NA_real_, cubature::adaptIntegrate(bvtnormcub(x,
      #   sigma_e_i = sqrt(sigma_hat_sq * M[i, i]),
      #   sigma_e_j = sqrt(sigma_hat_sq * M[j, j]),
      #   rho = M[i, j] / sqrt(M[i, i] * M[j, j])),
      #   lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))$integral -
      #   E0[i] * E0[j])))), nrow = n, ncol = n)

      Sigma <- matrix(data = unlist(lapply(1:n, function(i) lapply(1:n,
      function(j) ifelse(j <= i, NA_real_, cubature::adaptIntegrate(
        normexpect_integrand(xx, sigma1 = sqrt(sigma_hat_sq * M[i, i]),
        sigma2 = sqrt(sigma_hat_sq * M[j, j]), rho =  M[i, j] /
          sqrt(M[i, i] * M[j, j])),
      lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))$integral -
      E0[i] * E0[j])))), nrow = n, ncol = n)

      # Sigma <- matrix(data = unlist(lapply(1:n, function(i) lapply(1:n,
      # function(j) ifelse(j <= i, NA_real_, cubature::adaptIntegrate(
      # normexp_integrand2(x, Sigmat = sigma_hat_sq * M[c(i, j), c(i, j)]),
      # lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))$integral -
      # E0[i] * E0[j])))), nrow = n, ncol = n)

      if (any(is.infinite(Sigma))) stop("Numerical integration returned infinite values")

      Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
      diag(Sigma) <- sqrt(2 * sigma_hat_sq * diag(M)) / pi *
        (sqrt(pi) - gamma(3 / 4) ^ 2)
    }

    trfunc <- function(A, j) sum(diag(expm::`%^%`(A, j)))
    inside <- (Cmat - teststat * Bmat) %*% Sigma
    trvec <- vapply(1:3, function(j) trfunc(inside, j), NA_real_)
    teststatpval <- (trvec[2] ^ 3 - prod(trvec)) / trvec[3] ^ 2
    df <- trvec[2] ^ 3 / trvec[3] ^ 2
    pval <- stats::pchisq(teststatpval, df = df, lower.tail = FALSE)

  } else if (distmethod == "bootstrap") {

    Tstar <- rep(NA_real_, B)
    if (!is.null(seed)) set.seed(seed)
    for (b in 1:B) {
      ystar <- stats::rnorm(n, mean = yhat,
                     sd = sqrt(sigma_hat_sq))
      lmstar <- stats::lm.fit(X, ystar)
      estar <- lmstar$residuals
      sigma_hat_sqstar <- sum(estar ^ 2) / (n - p)
      E0star <- 2 ^ (1 / 4) * gamma(3 / 4) *
        sqrt(sigma_hat_sqstar * diag(M)) / sqrt(pi)
      sstar <- sqrt(abs(estar)) - E0star
      Tstar[b] <- as.double((t(sstar) %*% Cmat %*% sstar) /
        (t(sstar) %*% Bmat %*% sstar))
    }
    pval <- sum(Tstar >= teststat) / B
    df <- NULL

  } else stop("Invalid value of argument `method`")

  rval <- structure(list(statistic = teststat, p.value = pval, parameter = df,
               null.value = "Homoskedasticity", alternative = "greater"),
               class = "htest")
  broom::tidy(rval)
}
