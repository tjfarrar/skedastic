#' Probability mass function of nonparametric trend statistic \eqn{D}
#'
#' This function computes \eqn{\Pr(D = k)}, i.e. the probability mass
#' function for \eqn{D=\sum_{i=1}^{n} (R_i - i)^2}, the nonparametric trend
#' statistic proposed by \insertCite{Lehmann75;textual}{skedastic}, under the
#' assumption that the ranks \eqn{R_i} are computed on a series of \eqn{n}
#' independent and identically distributed random variables with no ties. The
#' function is used within \code{\link{horn}} in computing \eqn{p}-values for
#' Horn's nonparametric test for heteroskedasticity in a linear regression
#' model \insertCite{Horn81}{skedastic}. The support of \eqn{D}
#' consists of consecutive even numbers from 0 to \eqn{\frac{n(n-1)(n+1)}{3}},
#' with the exception of the case \eqn{n=3}, when the value 4 is excluded from
#' the support.
#'
#' @param n A positive integer representing the number of observations in the
#' series. Note that computation time increases rapidly with \eqn{n} and is
#' infeasible for \eqn{n>11}.
#' @param k An integer of \code{length} 1
#' or a character \code{"all"} indicating that the probability mass function
#' should be computed across the entire support of \eqn{D}.
#'
#' @return A \code{data.frame} containing two objects: \code{"value"}, an
#'   integer vector containing values specified in \code{k} (or the support of
#'   \eqn{D}, if \code{k == "all"}); and \code{"prob"}, a double vector
#'   representing the probability that \eqn{D} takes on each corresponding
#'   value in the \code{"value"} object
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{horn}}
#'
#' @examples
#' prob <- dDtrend(k = "all", n = 9)
#' values <- as.integer(names(dvar))
#' plot(c(values[1], values[1]), c(0, prob[1]), type = "l",
#'   axes = FALSE, xlab = expression(k), ylab = expression(Pr(D == k)),
#'   xlim = c(0, 250), yaxs = "i", ylim = c(0, 1.05 * max(prob)))
#'   axis(side = 1, at = seq(0, 250, 25), las = 2)
#' for (i in seq_along(values)) {
#'  lines(c(values[i], values[i]), c(0, prob[i]))
#' }
#'

dDtrend <- function(k = "all", n) {

  if (n > 10) warning("Computation of dDtrend is very slow for n > 10 and
                       computation time is of order n!")
  perms <- arrangements::permutations(n, n)
  valtab <- table(apply(perms, 1, function(x) sum((x - 1:n) ^ 2)))
  # value <- as.integer(names(valtab))
  prob <- as.double(valtab / factorial(n))
  names(prob) <- names(valtab)

  if (is.character(k)) {
    if (k != "all") stop("Invalid character value for k")
    return(prob)
  } else if (length(k) == 1) {
    if (k %% 1 != 0) stop("Invalid integer value for k")
    return(prob[as.integer(names(prob)) == k])
    # return(list("value" = k, "prob" = prob[value == k]))
  } else if (length(k) > 1) {
    stop("Argument `k` must be an integer of length 1 or a character, \"all\"")
  }
}

#' Cumulative distribution function of nonparametric trend statistic \eqn{D}
#'
#' This function computes \eqn{\Pr(D \le k)} (\eqn{\Pr(D \ge k)}), i.e.
#' lower (upper) cumulative probabilities for
#' \eqn{D=\sum_{i=1}^{n} (R_i - i)^2}, the nonparametric trend statistic
#' proposed by \insertCite{Lehmann75;textual}{skedastic}, under the assumption
#' that the ranks \eqn{R_i} are computed on a series of \eqn{n} independent and
#' identically distributed random variables with no ties. The function may be
#' used to compute one-sided \eqn{p}-values for the nonparametric test for
#' heteroskedasticity of \insertCite{Horn81;textual}{skedastic}. Computation
#' time is extremely slow for \eqn{n > 10} if \code{usedata} is set to
#' \code{FALSE}; thus \code{\link{horn}} uses a normal approximation in that
#' case.
#'
#' @param n A positive integer representing the number of observations in the
#'    series.
#' @param k An integer of \code{length} 1; the value must be between 0 and
#'    \eqn{\frac{n(n-1)(n+1)}{3}} inclusive.
#' @param lower.tail A logical. Should lower tailed cumulative probability be
#'    calculated? Defaults to \code{TRUE}. Note that both lower- and upper-
#'    tailed cumulative probabilities are computed inclusive of \code{k}.
#' @param exact A logical. Should exact distribution of \eqn{D} be used by
#'    calling \code{\link{dDtrend}}? If \code{FALSE}, a normal approximation
#'    is used. If \code{tiefreq} is not \code{NULL} (ties are present),
#'    normal approximation is used regardless of the value of \code{exact}.
#' @param correct A logical. Should continuity correction be used in normal
#'    approximation? Defaults to \code{TRUE}. Since the support of \eqn{D}
#'    consists always of even integers beginning from 0, the continuity
#'    correction consists of adding 1 to \eqn{D} when approximating the
#'    lower-tail probability \eqn{\Pr(D \le k)} and subtracting 1 from
#'    \eqn{D} when approximating the upper-tail probability
#'    \eqn{\Pr(D \ge k)}.
#' @param tiefreq A double vector corresponding to the value of \eqn{d_i}
#'    in \insertCite{Lehmann75;textual}{skedastic}. These are the frequencies
#'    of the various tied ranks. If ties are absent, \code{NULL} (the default).
#'
#' @return A double between 0 and 1 representing the probability/ies of \eqn{D}
#'    taking on at least (or at most) a value of \eqn{k}.
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{dDtrend}}, \code{\link{horn}}
#'
#' @examples
#' # For an independent sample of size 9, the probability that D is <= 50 is
#' # 0.05399857
#' pDtrend(k = 50, n = 9)
#' # Normal approximation of the above with continuity correction is
#' # 0.05193808
#' pDtrend(k = 50, n = 9, exact = FALSE, correct = TRUE)
#' # For an independent sample of size 50, the probability that D is >= 20000 is
#' # is 0.6093583
#' pDtrend(k = 2e4, n = 50, lower.tail = FALSE)
#'

pDtrend <- function(k, n, lower.tail = TRUE,
                    exact = (n <= 10), correct = TRUE, tiefreq = NULL) {

  if (exact && is.null(tiefreq)) {
    exact_upper <- function(j, m = n) {
      dvar <- dDtrend(n = m)
      values <- as.integer(names(dvar))
      return(sum(prob[values >= k]))
    }
    exact_lower <- function(j, m = n) {
      dvar <- dDtrend(n = m)
      values <- as.integer(names(dvar))
      return(sum(prob[values <= k]))
    }
    if (lower.tail) {
      vapply(k, exact_lower, FUN.VALUE = NA_real_)
    } else {
      vapply(k, exact_upper, FUN.VALUE = NA_real_)
    }
  } else {
    if (is.null(tiefreq)) {
      ED <- (n ^ 3 - n) / 6
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36
      normapprox_upper <- function(j, contcorrect = correct) {
        return(stats::pnorm((j - contcorrect - ED) / sqrt(VD), lower.tail = FALSE))
      }
      normapprox_lower <- function(j, contcorrect = correct) {
        return(stats::pnorm((j + contcorrect - ED) / sqrt(VD), lower.tail = TRUE))
      }
    } else {
      ED <- (n ^ 3 - n) / 6 - sum(tiefreq ^ 3 - tiefreq) / 12
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36 *
               (1 - sum(tiefreq ^ 3 - tiefreq) / (n ^ 3 - n))
      normapprox_upper <- function(j) {
        return(stats::pnorm((j - ED) / sqrt(VD), lower.tail = FALSE))
      }
      normapprox_lower <- function(j) {
        return(stats::pnorm((j - ED) / sqrt(VD), lower.tail = TRUE))
      }
    }
    if (lower.tail) {
      vapply(k, normapprox_lower, FUN.VALUE = NA_real_)
    } else {
      vapply(k, normapprox_upper, FUN.VALUE = NA_real_)
    }
  }
}
