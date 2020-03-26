#' Probability mass function of nonparametric trend statistic \eqn{D}
#'
#' This function computes \eqn{\Pr(D = k)}, i.e. the probability mass
#' function for \eqn{D=\sum_{i=1}^{n} (R_i - i)^2}, the nonparametric trend
#' statistic proposed by \insertCite{Lehmann75;textual}{skedastic}, under the
#' assumption that the ranks \eqn{R_i} are computed on a series of \eqn{n}
#' independent and identically distributed random variables with no ties. The
#' function is used within \code{\link{horn}} in computing \eqn{p}-values for
#' Horn's nonparametric test for heteroskedasticity in a linear regression model
#' (when sample size \eqn{n \leq 10}) \insertCite{Horn81}{skedastic}. Note that
#' Lehmann and Horn use a \eqn{D} statistic without the \eqn{\frac{1}{2}}
#' factor; this has been included so that support consists of whole numbers
#' incremented by 1, rather than incremented by 2. The support of \eqn{D}
#' consists of consecutive even numbers from 0 to \eqn{\frac{n(n-1)(n+1)}{3}},
#' with the exception of the case \eqn{n=3}, when the value 4 is excluded from
#' the support.
#'
#' @param n A positive integer representing the number of observations in the
#' series. Note that computation time increases rapidly with \eqn{n} and is
#' infeasible for \eqn{n>11}.
#' @param k An integer or a sequence of integers strictly incrementing by 1,
#'   or a character "all" indicating that the probability mass function should
#'   be computed across the entire support of \eqn{D}.
#'
#' @return A \code{data.frame} containing two objects: \code{"value"}, an
#'   integer vector containing values specified in \code{k} (or the support of
#'   \eqn{D}, if \code{k == "all"}); and \code{"prob"}, a double vector
#'   representing the probability that \eqn{D} takes on each corresponding
#'   value in the \code{"values"} object
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{horn}}
#'
#' @examples
#' dvar <- dDtrend(9)
#' plot(c(dvar$value[1], dvar$value[1]), c(0, dvar$prob[1]), type = "l",
#'      axes = FALSE, xlab = expression(k), ylab = expression(Pr(D == k)),
#'      xlim = c(0, 250), yaxs = "i", ylim = c(0, 1.05 * max(dvar$prob)))
#' axis(side = 1, at = seq(0, 250, 25), las = 2)
#' for (i in seq_along(dvar$value)) {
#'    lines(c(dvar$value[i], dvar$value[i]), c(0, dvar$prob[i]))
#' }
#'

dDtrend <- function(n, k = "all") {

  if (n > 10) warning("Computation of dDtrend is very slow for n > 10 and
                       computation time is of order n!")
  perms <- arrangements::permutations(n, n)
  valtab <- table(apply(perms, 1, function(x) sum((x - 1:n) ^ 2)))
  value <- as.integer(names(valtab))
  prob <- as.double(valtab / factorial(n))

  if (is.character(k)) {
    if (k != "all") stop("Invalid character value for k")
    return(list("value" = value, "prob" = prob))
  } else if (length(k) == 1) {
    if (k %% 1 != 0) stop("Invalid integer value for k")
    return(list("value" = k, "prob" = prob[value == k]))
  } else if (length(k) > 1) {
    if (!all(k[-1] - k[-length(k)]) == 1 || k[1] %% 1 != 0) stop("At least
                                          one invalid integer value for k")
    return(list("value" = k, "prob" = prob[value %in% k]))
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
#' @param k An integer or a sequence of integers strictly incrementing by 1,
#'    with all values between 0 and \code{n - 1} inclusive. Represents the
#'    number of peaks in the series.
#' @param lower.tail A logical. Should lower tailed cumulative probability be
#'    calculated? Defaults to \code{TRUE}. Note that both lower- and upper-
#'    tailed cumulative probabilities are computed inclusive of \code{k}.
#' @param exact A logical. Should exact distribution of \eqn{D} be used by
#'    calling \code{\link{dDtrend}}? If \code{FALSE}, a normal approximation
#'    is used.
#' @param correct A logical. Should continuity correction be used in normal
#'    approximation? Defaults to \code{TRUE}. Since the support of \eqn{D}
#'    consists always of even integers beginning from 0, the continuity
#'    correction consists of adding 1 to \eqn{D} when approximating the
#'    lower-tail probability \eqn{\Pr(D \le k)} and subtracting 1 from
#'    \eqn{D} when approximating the upper-tail probability
#'    \eqn{\Pr(D \ge k)}.
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
#' pDtrend(9, 50)
#' # Normal approximation of the above with continuity correction is
#' # 0.05193808
#' pDtrend(9, 50, exact = FALSE, correct = TRUE)
#' # For an independent sample of size 50, the probability that D is >= 20000 is
#' # is 0.6093583
#' pDtrend(50, 2e4, lower.tail = FALSE)
#'

pDtrend <- function(n, k, lower.tail = TRUE,
                    exact = (n <= 10), correct = !exact) {

  if (exact) {
    exact_upper <- function(j, m = n) {
      dvar <- dDtrend(m)
      return(sum(dvar$prob[dvar$value >= k]))
    }
    exact_lower <- function(j, m = n) {
      dvar <- dDtrend(m)
      return(sum(dvar$prob[dvar$value <= k]))
    }
    if (lower.tail) {
      vapply(k, exact_lower, FUN.VALUE = NA_real_)
    } else {
      vapply(k, exact_upper, FUN.VALUE = NA_real_)
    }
  } else {
    normapprox_upper <- function(j, m = n, contcorrect = correct) {
      ED <- (m ^ 3 - m) / 6
      VD <- (m ^ 2 * (m + 1) ^ 2 * (m - 1)) / 36
      return(stats::pnorm((j - contcorrect - ED) / sqrt(VD), lower.tail = FALSE))
    }
    normapprox_lower <- function(j, m = n, contcorrect = correct) {
      ED <- (m ^ 3 - m) / 6
      VD <- (m ^ 2 * (m + 1) ^ 2 * (m - 1)) / 36
      return(stats::pnorm((j + contcorrect - ED) / sqrt(VD), lower.tail = TRUE))
    }
    if (lower.tail) {
      vapply(k, normapprox_lower, FUN.VALUE = NA_real_)
    } else {
      vapply(k, normapprox_upper, FUN.VALUE = NA_real_)
    }
  }
}
