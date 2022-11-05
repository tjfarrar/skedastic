#' Probabilities for a Ratio of Quadratic Forms in a Normal Random Vector
#'
#' This function computes cumulative probabilities (lower or upper tail) on a
#'    ratio of quadratic forms in a vector of normally distributed
#'    random variables.
#'
#' @details Most of the work is done by other functions, namely
#'    \code{\link[CompQuadForm]{imhof}}, \code{\link[CompQuadForm]{davies}},
#'    or \code{\link[stats]{integrate}} (depending on the \code{algorithm}
#'    argument). It is assumed that the ratio of quadratic forms can be
#'    expressed as
#'    \deqn{R = \displaystyle\frac{x' A x}{x' B x}} where \eqn{x} is an
#'    \eqn{n}-dimensional normally distributed random variable with mean vector
#'    \eqn{\mu} and covariance matrix \eqn{\Sigma}, and \eqn{A} and
#'    \eqn{B} are real-valued, symmetric \eqn{n\times n} matrices. Matrix
#'    \eqn{B} must be non-negative definite to ensure that the denominator of
#'    the ratio of quadratic forms is nonzero.
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
#' @param r A double representing the value(s) for which \eqn{\Pr(R\le r)} or
#'    \eqn{\Pr(R \ge r)} should be computed.
#' @param A A numeric, symmetric matrix that is symmetric
#' @param B A numeric, symmetric, non-negative definite matrix having the same
#'    dimensions as \code{A}.
#' @param Sigma A numeric, symmetric matrix with the same dimensions as
#'    \code{A} and \code{B}, denoting the covariance matrix of the normal
#'    random vector. Defaults to the identity matrix, corresponding to the case
#'    in which the normal random variables are independent and identically
#'    distributed.
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
#' @param usenames A logical. If \code{TRUE}, the function value has a
#'    \code{names} attribute corresponding to \code{r}.
#'
#' @return A double denoting the probability/ies corresponding to the value(s)
#'    \code{r}.
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
#' pRQF(r = 1, A = A, B = B)
#' pRQF(r = 1, A = A, B = B, algorithm = "integrate")
#' pRQF(r = 1:3, A = A, B = B, algorithm = "davies")
#'

pRQF <- function(r, A, B, Sigma = diag(nrow(A)),
                    algorithm = c("imhof", "davies", "integrate"),
                    lower.tail = TRUE, usenames = FALSE) {

  algorithm <- match.arg(algorithm, c("imhof", "davies", "integrate"))

  n <- dim(A)
  if (!all.equal(n, dim(B)) || !all.equal(n, dim(Sigma)))
    stop("Dimension mismatch between A, B, and/or Sigma")
  if (!is.symmetric.mat(Sigma)) stop("Covariance matrix must be symmetric")
  if (all(B == 0))
    stop("B cannot be the zero matrix; this would make the denominator of the RQF zero")

   if (!is.pos.semidef.mat(B)) stop("B must be a positive semi-definite matrix")

  if (all(Sigma == diag(n[1]))) {
    Astar <- A
    Bstar <- B
  } else {
    Sigmasqrt <- pracma::sqrtm(Sigma)$B
    Astar <- Sigmasqrt %*% A %*% Sigmasqrt
    Bstar <- Sigmasqrt %*% B %*% Sigmasqrt
  }

  oner <- function(s) {
    PLambP <- Astar - s * Bstar

    mylambda <- bazar::almost.unique(eigen(PLambP, only.values = TRUE)$values)
    if (any(abs(Im(mylambda)) > sqrt(.Machine$double.eps)))
      stop("Astar - r * Bstar has complex eigenvalue(s)")
    mylambda <- mylambda[(abs(Re(mylambda)) > .Machine$double.eps &
                         abs(Im(mylambda)) < .Machine$double.eps)]
    mylambda <- as.numeric(Re(mylambda))

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
  vapply(r, oner, NA_real_, USE.NAMES = usenames)
}
