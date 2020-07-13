#' Computation of Conditional Two-Sided \eqn{p}-Values
#'
#' Computes the conditional \eqn{p}-value \eqn{P_C} for a continuous
#'    or discrete test statistic, as defined in
#'    \insertCite{Kulinskaya08;textual}{skedastic}. This provides a method
#'    for computing a two-sided \eqn{p}-value from an asymmetric null
#'    distribution.
#'
#' Let \eqn{T} be a statistic that, under the null hypothesis, has cumulative
#'    distribution function \eqn{F} and probability density or mass function
#'    \eqn{f}. Denote by \eqn{A} a generic location parameter chosen to separate
#'    the two tails of the distribution. Particular examples include the mean
#'    \eqn{E(T|\mathrm{H}_0)}, the mode \eqn{\arg \sup_{t} f(t)}, or the median
#'    \eqn{F^{-1}\left(\frac{1}{2}\right)}. Let \eqn{q} be the observed value
#'    of \eqn{T}.
#'
#'    In the continuous case, the conditional two-sided \eqn{p}-value centered
#'    at \eqn{A} is defined as
#'    \deqn{P_C^A(q)=\frac{F(q)}{F(A)}1_{q \le A} + \frac{1-F(q)}{1-F(A)}1_{q > A}}
#'    where \eqn{1_{\cdot}} is the indicator function. In the discrete case,
#'    \eqn{P_C^A} depends on whether \eqn{A} is an attainable value within the
#'    support of \eqn{T}. If \eqn{A} is not attainable, the conditional two-sided
#'    \eqn{p}-value centred at \eqn{A} is defined as
#'    \deqn{P_C^{A}(q)=\frac{\Pr(T\le q)}{\Pr(T<A)}1_{q<A} + \frac{\Pr(T\ge q)}{\Pr(T>A)}1_{q>A}}
#'    If \eqn{A} is attainable, the conditional two-sided \eqn{p}-value centred
#'    at \eqn{A} is defined as
#'    \deqn{P_C^{A}(q)=\frac{\Pr(T\le q)}{\Pr(T\le A)/\left(1+\Pr(T=A)\right)} 1_{q<A} +
#'    1_{q=A}+\frac{\Pr(T\ge q)}{\Pr(T \ge A)/\left(1+\Pr(T=A)\right)} 1_{q>A}}
#'
#' @param q A double representing the quantile, i.e. the observed value of the
#'    test statistic for which a two-sided \eqn{p}-value is to be computed
#' @param Aloc a double representing a generic location parameter chosen to
#'    separate the tails of the distribution. Note that if \code{Aloc}
#'    corresponds to the median of \code{CDF}, there is no difference between
#'    the two methods, in the continuous case. If \code{Aloc} is not specified,
#'    the function attempts to compute the expectation of \code{CDF} using
#'    numerical integration and, if successful, uses this as \code{Aloc}.
#'    However, this may yield unexpected results, especially if \code{CDF} is
#'    not one of the cumulative distribution functions of well-known
#'    distributions included in the \code{stats} package.
#' @param CDF A function representing the cumulative distribution function of
#'    the test statistic under the null hypothesis, i.e.
#'    \eqn{\Pr(T\le q|\mathrm{H}_0)}.
#' @param method A character specifying the method to use to calculate
#'    two-sided \eqn{p}-value; one of \code{"doubled"} (representing
#'    doubling of the one-sided \eqn{p}-value),
#'    \code{"kulinskaya"} (representing the method of
#'    \insertCite{Kulinskaya08;textual}{skedastic}), or \code{"minlikelihood"}
#'    (representing the sum of probabilities for values with probability less
#'    than or equal to that of the observed value. Partial matching is used.
#'    Note that the \code{"minlikelihood"} method is available only
#'    for discrete distributions.
#' @param continuous A logical indicating whether the test statistic is a
#'    continuous (\code{TRUE}) or discrete (\code{FALSE}) random variable.
#'    Defaults to \code{TRUE}.
#' @param supportlim A numeric vector of \code{length} 2, giving the minimum
#'    and maximum values in the support of the distribution whose cumulative
#'    distribution function is \code{CDF}. This argument is only used if the
#'    distribution is discrete (i.e. if \code{continuous} is \code{FALSE}) and
#'    if \code{method} is \code{"minlikelihood"} or \code{Aloc} is not
#'    specified. If \code{supportlim} is not supplied, the function assumes
#'    that the support is \code{-1e6:1e6}. Values of \code{-Inf} and \code{Inf}
#'    may be supplied, but if so, the support is truncated at \code{-1e6} and
#'    \code{1e6} respectively.
#' @param ... Optional arguments to pass to \code{CDF}.
#'
#' @return A double.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' # Computation of two-sided p-value for F test for equality of variances
#' n1 <- 10
#' n2 <- 20
#' set.seed(1234)
#' x1 <- stats::rnorm(n1, mean = 0, sd = 1)
#' x2 <- stats::rnorm(n2, mean = 0, sd = 3)
#' # 'Conventional' two-sided p-value obtained by doubling one-sided p-value:
#' stats::var.test(x1, x2, alternative = "two.sided")$p.value
#' # This is replicated in `twosidedpval` by setting `method` argument to `"doubled"`
#' twosidedpval(q = var(x1) / var(x2), CDF = stats::pf, continuous = TRUE,
#'  method = "doubled", Aloc = 1, df1 = n1 - 1, df2 = n2 - 1)
#' # Conditional two-sided p-value centered at df (mean of chi-squared r.v.):
#' twosidedpval(q = var(x1) / var(x2), CDF = stats::pf, continuous = TRUE,
#'  method = "kulinskaya", Aloc = 1, df1 = n1 - 1, df2 = n2 - 1)

twosidedpval <- function(q, CDF, continuous, method = c("doubled", "kulinskaya",
                         "minlikelihood"), Aloc, supportlim, ...) {

  method <- match.arg(method, c("doubled", "kulinskaya", "minlikelihood"))

  if (missing(Aloc) && method != "minlikelihood") {
    Aloc <- meanfromCDF(theCDF = CDF, cont = continuous, suplim = supportlim,
                        ...)
  }

  if (method == "kulinskaya") {
    if (continuous) {
      CDF(q, ...) / CDF(Aloc, ...) * (q <= Aloc) + (1 - CDF(q, ...)) /
        (1 - CDF(Aloc, ...)) * (q >= Aloc)
    } else {
      if (!(value_possible(x = Aloc, myCDF = CDF, ...))) {
        wL <- CDF(Aloc, ...)
        wR <- 1 - CDF(Aloc, ...)
        CDF(q, ...) / wL * (q < Aloc) + (q == Aloc) + (1 - CDF(q - 1, ...)) /
          wR * (q > Aloc)
      } else {
        wLm <- CDF(Aloc, ...) / (1 + CDF(Aloc, ...) - CDF(Aloc - 1, ...))
        wRm <- (1 - CDF(Aloc - 1, ...)) / (1 + CDF(Aloc, ...) - CDF(Aloc - 1, ...))
        CDF(q, ...) / wLm * (q < Aloc) + (q == Aloc) + (1 - CDF(q - 1, ...)) /
          wRm * (q > Aloc)
      }
    }
  } else if (method == "doubled") {
    min(2 * CDF(q, ...), 2 * (1 - CDF(q, ...)))
  } else if (method == "minlikelihood") {
    if (continuous) {
      stop("Method minlikelihood not implemented for continuous function")
    } else {
      if (missing(supportlim)) {
        support <- -1e6:1e6
      } else {
        support <- supportlim[1]:supportlim[2]
      }
      allPMF <- vapply(support, function(j) CDF(j, ...) -
                         CDF(j - 1, ...), NA_real_)
      PMFq <- CDF(q, ...) - CDF(q - 1, ...)
      sum(allPMF[allPMF <= PMFq])
    }
  } else stop("Invalid `method` argument")
}
