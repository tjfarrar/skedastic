#' Count peaks in a data series
#'
#' This function computes the number of peaks in a univariate data series, with
#' peak defined as per \insertCite{Goldfeld65;textual}{skedastic}. The function
#' is used in the Goldfeld-Quandt nonparametric test for heteroskedasticity in
#' a linear model.
#'
#' @param x A double vector.
#'
#' @return An integer value between \code{0} and \code{length(x) - 1}
#'    representing the number of peaks in the series.
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{goldfeld_quandt}}
#'
#' @examples
#' set.seed(9586)
#' countpeaks(rnorm(20))
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' countpeaks(mtcars_lm$residuals)
#'

countpeaks <- function(x) {
  n <- length(x)
  npeaks <- 0L
  peak <- x[1]
  for (i in 2:n) {
    if (x[i] >= peak) {
      npeaks <- npeaks + 1L
      peak <- x[i]
    }
  }
  return(npeaks)
}

#' Probability mass function of number of peaks in a random series
#'
#' This function computes \eqn{P(n,k)} as defined by
#' \insertCite{Goldfeld65;textual}{skedastic}, i.e. the probability that a
#' series of \eqn{n} independent and identically distributed random variables
#' contains exactly \eqn{k} peaks, with peaks also as defined by
#' \insertCite{Goldfeld65;textual}{skedastic}. The function is used in
#' \code{\link{ppeak}} to compute \eqn{p}-values for the Goldfeld-Quandt
#' nonparametric test for heteroskedasticity in a linear model.
#'
#' @param n A positive integer representing the number of observations in the
#' series.
#' @param k An integer or a sequence of integers strictly incrementing by 1,
#' with all values between 0 and \code{n - 1} inclusive. Represents the number
#' of peaks in the series.
#' @param usedata A logical. Should probability mass function values be
#' read from \code{\link{dpeakdat}} rather than computing them? This option
#' will save significantly on computation time if \eqn{n < 170} but is
#' currently only available for \eqn{n \leq 1000}.
#'
#' @return A double between 0 and 1 representing the probability of exactly
#' k peaks occurring in a series of \eqn{n} independent and identically
#' distributed continuous random variables. Computation time is very slow for
#' \eqn{n > 170} (if \code{usedata} is \code{FALSE}) and for \eqn{n > 1000}
#' regardless of \code{usedata} value.
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{ppeak}}, \code{\link{goldfeld_quandt}}
#'
#' @examples
#' dpeak(10, 0:9)
#' plot(0:9, dpeak(10, 0:9), type = "p", pch = 20, xlab = "Number of Peaks",
#'          ylab = "Probability")
#'
#' # `dpeakdat` is a dataset containing probabilities generated from `dpeak`
#' utils::data(dpeakdat)
#' expval <- unlist(lapply(dpeakdat,
#'                  function(p) sum(p * 0:(length(p) - 1))))
#' plot(1:1000, expval[1:1000], type = "l", xlab = parse(text = "n"),
#'      ylab = "Expected Number of Peaks")

dpeak <- function(n, k, usedata = FALSE) {
  if (length(k) == 1) {
    maxk <- k
  } else if (length(k) > 1) {
    if (!all(k[-1] - k[-length(k)]) == 1) stop("Only sequential k values permitted")
    maxk <- max(k)
  }
  if (maxk >= n) {
    warning(paste0("No. of peaks `k` must be less than no. of observations `n`\n
            k truncated at n - 1 = ", n - 1))
    k <- ifelse(length(k) == 1, n - 1, min(k):(n - 1))
    maxk <- n - 1
  }
  if (usedata) {
    utils::data(dpeakdat)
    dpeakdat[[n]][k+1]
  } else {
    if (n > 170) {
      N <- Rmpfr::mpfrArray(0, precBits = 53, dim = c(n, maxk + 1))
      factorial_denom <- gmp::factorialZ(n)
    } else {
      N <- matrix(data = 0, nrow = n, ncol = maxk + 1) # k can be 0
      factorial_denom <- factorial(n)
    }
    N[1, 1] <- 1
    if (n > 1) {
      for (i in 2:n) {
        N[i, 1] <- (i - 1) * N[i - 1, 1]
        if (maxk >= 1) {
          for (j in 2:(maxk + 1)) {
            N[i, j] <- (i - 1) * N[i - 1, j] + N[i - 1, j - 1]
          }
        }
      }
    }
    as.double(N[n, k + 1] / factorial_denom)
  }
}

#' Cumulative distribution function of number of peaks in a random series
#'
#' This function computes \eqn{\sum_{k} P(n,k)}, i.e. the probability that a
#' series of \eqn{n} independent and identically distributed random variables
#' contains \eqn{\ge k} \eqn{(\le k)} peaks, with peaks as defined in
#' \insertCite{Goldfeld65;textual}{skedastic}. The function may be used to
#' compute \eqn{p}-values for the Goldfeld-Quandt nonparametric test for
#' heteroskedasticity in a linear model. Computation time is very slow for
#' \eqn{n > 170} if \code{usedata} is set to \code{FALSE}.
#'
#' @param n A positive integer representing the number of observations in the
#'    series.
#' @param k An integer or a sequence of integers strictly incrementing by 1,
#'    with all values between 0 and \code{n - 1} inclusive. Represents the
#'    number of peaks in the series.
#' @param upper A logical. Should upper tailed cumulative probability be
#'    calculated? Defaults to \code{TRUE} due to function being designed
#'    primarily for calculating \eqn{p}-values for the peaks test. Note that
#'    both upper and lower tailed cumulative probabilities are computed
#'    inclusive of \code{k}.
#' @param usedata A logical. Should probability mass function values be
#'    read from \code{\link{dpeakdat}} rather than computing them from
#'    \code{\link{dpeak}}? This option will save significantly on
#'    computation time if \eqn{n < 170} but is currently only available
#'    for \eqn{n \le 1000}.
#'
#' @return A double between 0 and 1 representing the probability of at least
#'    (at most) k peaks occurring in a series of \eqn{n} independent and
#'    identically distributed continuous random variables.
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{dpeak}}, \code{\link{goldfeld_quandt}}
#'
#' @examples
#' # For an independent sample of size 250, the probability of at least 10
#' # peaks is 0.06186582
#' ppeak(250, 10, upper = TRUE, usedata = TRUE)
#' # For an independent sample of size 10, the probability of at most 2 peaks
#' # is 0.7060615
#' ppeak(10, 2, upper = FALSE, usedata = FALSE)
#'

ppeak <- function(n, k, upper = TRUE, usedata = TRUE) {

  if (usedata) {
    utils::data(dpeakdat)
    maxn_indata <- length(dpeakdat)
  }
  ppeak_up <- function(n, j) {
    if (j >= n) {
      stop("`k` must be less than `n`")
    } else if (j == 0) {
      return(1)
    } else {
      if (usedata && n <= maxn_indata) {
        return(sum(dpeakdat[[n]][(j + 1):(n)]))
      } else {
        return(sum(dpeak(n, j:(n - 1))))
      }
    }
  }

  ppeak_lo <- function(n, j) {
    if (j >= n) {
      stop("`k` must be less than `n`")
    } else {
      if (usedata && n <= maxn_indata) {
        return(sum(dpeakdat[[n]][1:(j + 1)]))
      } else {
        return(sum(dpeak(n, 0:j)))
      }
    }
  }

  if (upper) {
    vapply(k, ppeak_up, FUN.VALUE = NA_real_, n = n)
  } else {
    vapply(k, ppeak_lo, FUN.VALUE = NA_real_, n = n)
  }
}
