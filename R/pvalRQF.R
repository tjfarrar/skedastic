#' Probabilities for a Ratio of Quadratic Forms in i.i.d. Normal Random Variables
#'
#' This function computes cumulative probabilities (lower or upper tail) on a
#'    ratio of quadratic forms in a vector of i.i.d. normally distributed
#'    random variables. Most of the work is done by other functions, namely
#'    \code{\link[CompQuadForm]{imhof}}, \code{\link[CompQuadForm]{davies}},
#'    or \code{\link[stats]{integrate}} (depending on the \code{algorithm}
#'    argument). It is assumed that the ratio of quadratic forms can be
#'    expressed as
#'    \deqn{R = \displaystyle\frac{x' A x}{x' B x}} where \eqn{x} is an
#'    \eqn{n}-dimensional normally distributed random variable with mean vector
#'    0 and covariance matrix proportional to \eqn{I_n}, and \eqn{A} and
#'    \eqn{B} are real-valued, symmetric \eqn{n\times n} matrices.
#'
#' The function makes use of the fact that a probability statement involving a
#'    ratio of quadratic forms can be rewritten as a probability statement
#'    involving a quadratic form. Hence, methods for computing probabilities
#'    for a quadratic form in normal random variables, such as the Imhof
#'    algorithm \insertCite{Imhof61}{skedastic} or the Davies algorithm
#'    \insertCite{Davies80}{skedastic} can be applied to the rearranged
#'    expression to obtain the probability for the ratio of quadratic forms.
#'    Note that the Ruben-Farebrother algorithm (as implemented in
#'    \code{\link[CompQuadForm]{farebrother}}) cannot be used here because the
#'    \eqn{A} matrix within the quadratic form (after rearrangement of the
#'    probability statement involving a ratio of quadratic forms) is not in
#'    general positive semi-definite.
#'
#' The function throws a \code{warning} if the matrix \code{A - r * B} is not
#'    numerically symmetric, as determined by
#'    \code{\link[base]{isSymmetric.matrix}}. In such cases, the user should
#'    verify whether this is merely an issue of numerical precision or whether
#'    the matrix is actually not symmetric.
#'
#' @param r A double representing the value(s) for which \eqn{\Pr(R\le r)} or
#'    \eqn{\Pr(R \ge r)} should be computed.
#' @param A A numeric matrix that is symmetric
#' @param B A numeric matrix that is symmetric
#' @param algorithm A character, either \code{"imhof"}, \code{"davies"}, or
#'    \code{"integrate"}. Values \code{"imhof"} and \code{"integrate"}
#'    both implement the Imhof algorithm. The difference is that \code{"imhof"}
#'    means that \code{\link[CompQuadForm]{imhof}} is used, whereas
#'    \code{"integrate"} means that \code{\link[stats]{integrate}} is
#'    used (which is slower). The Imhof algorithm is more precise than the
#'    Davies algorithm.
#' @param lower.tail A logical. If \code{TRUE}, the cumulative distribution
#'    function \eqn{\Pr(R \le r)} is computed; if \code{FALSE}, the survival
#'    function \eqn{\Pr(R \ge r)} is computed.
#'
#' @return A double.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \insertCite{Duchesne10;textual}{skedastic}, the article associated
#'    with the \code{\link[CompQuadForm]{imhof}} and
#'    \code{\link[CompQuadForm]{davies}} functions.
#'
#' @examples
#' n <- 20
#' A <- matrix(data = 1, nrow = n, ncol = n)
#' B <- diag(n)
#' pvalRQF(r = 1, A = A, B = B)
#' pvalRQF(r = 1, A = A, B = B, algorithm = "integrate")
#' pvalRQF(r = 1:3, A = A, B = B, algorithm = "davies")
#'

pvalRQF <- function(r, A, B,
                    algorithm = c("imhof", "davies", "integrate"),
                    lower.tail = TRUE) {

  algorithm <- match.arg(algorithm, c("imhof", "davies", "integrate"))

  oner <- function(s) {
    PLambP <- A - s * B
    if(!isSymmetric.matrix(PLambP)) warning("Matrix A - r * B is not numerically symmetric")

    mylambda <- eigen(PLambP, only.values = TRUE)$values

    if (algorithm %in% c("imhof", "davies")) {

      upperprob <- do.call(`::`, args = list("CompQuadForm", algorithm))(q = 0,
                            lambda = mylambda)$Qq
      ifelse(lower.tail, 1 - upperprob, upperprob)
    } else if (algorithm == "integrate") {
      integrand <- function(u) {
        if (length(u) == 1) {
          theta <- 1 / 2 * sum(atan(mylambda * u))
          rho <- prod((1 + mylambda ^ 2 * u ^ 2) ^ (1 / 4))
          return(sin(theta) / (u * rho))
        } else if (length(u) > 1) {
          theta <- vapply(u, function(v) 1 / 2 * sum(atan(mylambda * v)), NA_real_)
          rho <- vapply(u, function(v) prod((1 + mylambda ^ 2 * v ^ 2) ^ (1 / 4)),
                        NA_real_)
          return(sin(theta) / (u * rho))
        }
      }
      upperprob <- 1 / pi * do.call(`::`, args = list("stats",
                  "integrate"))(f = integrand, lower = 0, upper = Inf)$value
      ifelse(lower.tail, 1 / 2 - upperprob, 1 / 2 + upperprob)
    }
  }
  vapply(r, oner, NA_real_, USE.NAMES = FALSE)
}
