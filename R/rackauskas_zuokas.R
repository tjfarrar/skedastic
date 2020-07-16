#' Rackauskas-Zuokas Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods of
#'    \insertCite{Rackauskas07;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' Rackauskas and Zuokas propose a class of tests that entails determining the
#'    largest weighted difference in variance of estimated error. The
#'    asymptotic behaviour of their test statistic \eqn{T_{n,\alpha}} is
#'    studied using the empirical polygonal process constructed from partial
#'    sums of the squared residuals. The test is right-tailed.
#' @param alpha A double such that \eqn{0 \le \alpha < 1/2}; a hyperparameter
#'    of the test. Defaults to 0.
#' @param pvalmethod A character, either \code{"data"} or \code{"sim"},
#'    determining which method to use to compute the empirical
#'    \eqn{p}-value. If \code{"data"}, the dataset \code{\link{T_alpha}}
#'    consisting of pre-generated Monte Carlo replicates from the
#'    asymptotic null distribution of the test statistic is loaded and used to
#'    compute empirical \eqn{p}-value. This is only available for certain
#'    values of \code{alpha}, namely \eqn{i/32} where \eqn{i=0,1,\ldots,15}.
#'    If \code{"sim"}, Monte Carlo replicates are generated from the
#'    asymptotic null distribution. Partial matching is used.
#' @param R An integer representing the number of Monte Carlo replicates to
#'    generate, if \code{pvalmethod == "sim"}. Ignored if
#'    \code{pvalmethod == "data"}.
#' @param m An integer representing the number of standard normal variates to
#'    use when generating the Brownian Bridge for each replicate, if
#'    \code{pvalmethod == "sim"}. Ignored if \code{pvalmethod == "data"}. If
#'    number of observations is small,
#'    \insertCite{Rackauskas07;textual}{skedastic} recommends using \eqn{m=n}.
#'    The dataset \code{\link{T_alpha}} used \eqn{m=2^17} which is
#'    computationally intensive.
#' @param sqZ A logical. If \code{TRUE}, the standard normal variates used
#'    in the Brownian Bridge when generating from the asymptotic null
#'    distribution are first squared, i.e. transformed to \eqn{\chi^2(1)}
#'    variates. This is recommended by
#'    \insertCite{Rackauskas07;textual}{skedastic} when the number of
#'    observations is small. Ignored if \code{pvalmethod == "data"}.
#' @param seed An integer representing the seed to be used for pseudorandom
#'    number generation when simulating values from the asymptotic null
#'    distribution. This is to provide reproducibility of test results.
#'    Ignored if \code{pvalmethod == "data"}. If user does not wish to set
#'    the seed, pass \code{NULL}.
#'
#' @inheritParams breusch_pagan
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
#' rackauskas_zuokas(mtcars_lm)
#' rackauskas_zuokas(mtcars_lm, alpha = 7 / 16)
#' \donttest{
#' n <- length(mtcars_lm$residuals)
#' rackauskas_zuokas(mtcars_lm, pvalmethod = "sim", m = n, sqZ = TRUE)
#' }
#'

rackauskas_zuokas <- function(mainlm, alpha = 0, pvalmethod = c("data", "sim"),
                               R = 2 ^ 14, m = 2 ^ 17, sqZ = FALSE, seed = 1234,
                              statonly = FALSE) {

  if (alpha < 0 || alpha >= 1 / 2) stop("Invalid `alpha` argument. `alpha` must be >= 0
              and < 1/2")

  processmainlm(m = mainlm, needy = FALSE, needp = FALSE)

  n <- length(e)
  Tnalpha <- max(vapply(1:(n - 1), function(ell) max((ell / n) ^ (-alpha) *
                      vapply(0:(n - ell), function(k)
      abs(sum(e[(k + 1):(k + ell)] ^ 2 - 1 / n * sum(e ^ 2))), NA_real_)), NA_real_))
  deltahat <- mean((e ^ 2 - mean(e ^ 2)) ^ 2)
  teststat <- Tnalpha / sqrt(deltahat * n)
  if (statonly) return(teststat)

  pvalmethod <- match.arg(pvalmethod, c("data", "sim"))
  if (pvalmethod  == "data") {
    utils::data(T_alpha)
    if (min(abs(alpha - (0:15 / 32))) > 1e-6) stop("Values from the
       null distribution have not been pre-generated for this
                                  value of alpha")
    whichcol <- alpha * 32 + 1
    pval <- sum(teststat < T_alpha[, whichcol]) / nrow(T_alpha)
  } else if (pvalmethod == "sim") {
    Talphavals <- rksim(R. = R, m. = m, sqZ. = sqZ, seed. = seed,
                        alpha. = alpha)
    pval <- sum(teststat < Talphavals) / R
  }
  rval <- structure(list(statistic = teststat, p.value = pval,
               null.value = "Homoskedasticity", alternative = "greater",
               method = pvalmethod, parameter = alpha), class = "htest")
  broom::tidy(rval)
}
