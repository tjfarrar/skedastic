#' Probability mass function of nonparametric trend statistic \eqn{D}
#'
#' This function computes \eqn{\Pr(D = k)}, i.e. the probability mass
#'   function for \eqn{D=\sum_{i=1}^{n} (R_i - i)^2}, the nonparametric trend
#'   statistic proposed by \insertCite{Lehmann75;textual}{skedastic}, under the
#'   assumption that the ranks \eqn{R_i} are computed on a series of \eqn{n}
#'   independent and identically distributed random variables with no ties.
#'
#' @details The function is used within \code{\link{horn}} in computing
#'   \eqn{p}-values for Horn's nonparametric test for heteroskedasticity in a
#'   linear regression model \insertCite{Horn81}{skedastic}. The support of
#'   \eqn{D} consists of consecutive even numbers from 0 to
#'   \eqn{\frac{n(n-1)(n+1)}{3}}, with the exception of the case \eqn{n=3},
#'   when the value 4 is excluded from the support. Note that computation speed
#'   for \code{k = "all"} is about the same as when \code{k} is set to an
#'   individual integer value, because the entire distribution is still
#'   computed in the latter case.
#'
#' @param n A positive integer representing the number of observations in the
#'   series. Note that computation time increases rapidly with \eqn{n} and is
#'   infeasible for \eqn{n>11}.
#' @param k An integer of \code{length} \eqn{\ge 1}
#'   or a character \code{"all"} (the default) indicating that the probability
#'   mass function should be applied to the entire support of \eqn{D}.
#' @param override A logical. By default, the function aborts if \eqn{n > 11}
#'   due to the prohibitively slow computation (which may cause some systems
#'   to crash). Setting this argument to \code{TRUE} overrides the abort.
#'
#' @return A double vector containing the probabilities corresponding to the
#'   integers in its \code{names} attribute.
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso \code{\link{horn}}
#'
#' @examples
#' prob <- dDtrend(k = "all", n = 9)
#' values <- as.integer(names(prob))
#' plot(c(values[1], values[1]), c(0, prob[1]), type = "l",
#'   axes = FALSE, xlab = expression(k), ylab = expression(Pr(D == k)),
#'   xlim = c(0, 250), yaxs = "i", ylim = c(0, 1.05 * max(prob)))
#'   axis(side = 1, at = seq(0, 250, 25), las = 2)
#' for (i in seq_along(values)) {
#'  lines(c(values[i], values[i]), c(0, prob[i]))
#' }
#'

dDtrend <- function(k = "all", n, override = FALSE) {

  kall <- FALSE
  if (is.na(k[1]) || is.null(k)) {
    stop("Argument k cannot be NA or NULL")
  } else if (is.character(k)) {
    if (length(k) == 1 && k == "all") {
      kall <- TRUE
    } else {
      stop("Invalid character value for k")
    }
  } else if (any(k %% 1 != 0)) {
    stop("Invalid value(s) in k; try an integer or \"all\"")
  }

  if (is.na(n[1]) || is.null(n)) {
    stop("Argument n cannot be NA or NULL")
  } else if (!is.numeric(n)) {
    stop("Argument n must be numeric")
  } else if (n %% 1 != 0) {
    stop("Invalid value for n; try an integer")
  } else if (n > 11 && !override) {
    stop("Computation of dDtrend is prohibitively slow for n > 11. Operation aborted. If user insists on proceeding, call function again with `override` set to `TRUE`.")
  }

  if (!requireNamespace("arrangements", quietly = TRUE)) {
    stop(
      "Package \"arrangements\" must be installed to use this function.",
      call. = FALSE
    )
  }

  perms <- arrangements::permutations(n, n)
  valtab <- table(apply(perms, 1, function(x) sum((x - 1:n) ^ 2)))
  prob <- as.double(valtab / factorial(n))
  names(prob) <- names(valtab)
  support <- as.integer(names(prob))
  if (kall) {
    return(prob)
  } else {
    if (any(!(k %in% support))) warning("One or more values in k not part of support of distribution")
    return(prob[support %in% k])
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
#' \code{FALSE}; thus a normal approximation is implemented, including a
#' continuity correction.
#'
#' @param n A positive integer representing the number of observations in the
#'    series.
#' @param k An integer of \code{length} \eqn{\ge 1}
#'   or a character \code{"all"} indicating that the cumulative
#'   distribution function should be applied to the entire support of \eqn{D}.
#'   The latter is only acceptable when \code{exact} is \code{TRUE}, since
#'   the distribution is otherwise continuous.
#' @param lower.tail A logical. Should lower tailed cumulative probability be
#'    calculated? Defaults to \code{TRUE}. Note that both lower- and upper-
#'    tailed cumulative probabilities are computed inclusive of \code{k}.
#' @param exact A logical. Should exact distribution of \eqn{D} be used by
#'    calling \code{\link{dDtrend}}? If \code{FALSE}, a normal approximation
#'    is used. If \code{tiefreq} is not \code{NA} (ties are present),
#'    normal approximation is used regardless of the value of \code{exact}.
#'    By default, \code{exact} is set to \code{TRUE} provided that
#'    \code{n <= 10}. Setting \code{exact} to \code{TRUE} for \code{n > 11}
#'    results in an error unless \code{override} is set to \code{TRUE}.
#' @param tiefreq A double vector corresponding to the value of \eqn{d_i}
#'    in \insertCite{Lehmann75;textual}{skedastic}. These are the frequencies
#'    of the various tied ranks. If ties are absent, \code{NA} (the default).
#' @param override A logical. By default, the \code{dDtrend} function aborts if
#'   \eqn{n > 11} due to the prohibitively slow computation (which may cause
#'   some systems to crash). Setting this argument to \code{TRUE} overrides
#'   the abort. Ignored unless \code{exact} is \code{TRUE}.
#'
#' @return A double between 0 and 1 representing the probability/ies of \eqn{D}
#'    taking on at least (at most) the value(s) in the \code{names} attribute.
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
#' pDtrend(k = 50, n = 9, exact = FALSE)
#' # For an independent sample of size 50, the probability that D is >= 20000 is
#' # is 0.6093583
#' pDtrend(k = 2e4, n = 50, lower.tail = FALSE)
#'

pDtrend <- function(k, n, lower.tail = TRUE, exact = (n <= 10), tiefreq = NA,
                    override = FALSE) {

  kall <- FALSE
  if (any(is.na(k)) || is.null(k)) {
    stop("Invalid k value(s)")
  } else if (is.character(k)) {
    if (length(k) == 1 && k == "all") {
      kall <- TRUE
    } else {
      stop("Invalid k value(s)")
    }
  } else if (any(!is.numeric(k))) stop("Invalid k value(s)")

  if (exact && (is.na(tiefreq[1]) || is.null(tiefreq))) {

    prob <- dDtrend(k = "all", n = n, override = override)
    values <- as.integer(names(prob))
    ineqfunc <- ifelse(lower.tail, `<=`, `>=`)
    if (kall) {
      cumprob <- vapply(values, function(j) sum(prob[ineqfunc(values, j)]),
                        NA_real_)
    } else {
      cumprob <- vapply(k, function(j) sum(prob[ineqfunc(values, j)]),
                        NA_real_)
    }

  } else {

    if (kall) stop("k = \"all\" is only valid when the exact, discrete distribution is used")
    if (is.na(tiefreq[1]) || is.null(tiefreq)) {
      ED <- (n ^ 3 - n) / 6
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36
      cumprob <- vapply(k, function(j) stats::pnorm((j - 1 - ED) / sqrt(VD),
                  lower.tail = lower.tail), NA_real_)
    } else {
      ED <- (n ^ 3 - n) / 6 - sum(tiefreq ^ 3 - tiefreq) / 12
      VD <- (n ^ 2 * (n + 1) ^ 2 * (n - 1)) / 36 *
               (1 - sum(tiefreq ^ 3 - tiefreq) / (n ^ 3 - n))
      cumprob <- vapply(k, function(j) stats::pnorm((j - ED) / sqrt(VD),
                  lower.tail = lower.tail), NA_real_)
    }
  }
  if (kall) {
    names(cumprob) <- names(prob)
  } else {
    names(cumprob) <- as.character(k)
  }
  cumprob
}
