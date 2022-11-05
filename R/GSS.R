#' Golden Section Search for Minimising Univariate Function over a Closed
#'    Interval
#'
#' Golden Section Search (GSS) is a useful algorithm for minimising a
#'    continuous univariate function \eqn{f(x)} over an interval
#'    \eqn{\left[a,b\right]} in instances where the first derivative
#'    \eqn{f'(x)} is not easily obtainable, such as with loss functions
#'    that need to be minimised under cross-validation to tune a
#'    hyperparameter in a machine learning model. The method is described
#'    by \insertCite{Fox21;textual}{skedastic}.
#'
#' @details This function is modelled after a MATLAB function written by
#'    \insertCite{Zarnowiec22;textual}{skedastic}. The method assumes that
#'    there is one local minimum in this interval. The solution produced is
#'    also an interval, the width of which can be made arbitrarily small by
#'    setting a tolerance value. Since the desired solution is a single
#'    point, the midpoint of the final interval can be taken as the best
#'    approximation to the solution.
#'
#'    The premise of the method is to shorten the interval
#'    \eqn{\left[a,b\right]} by comparing the values of the function at two
#'    test points, \eqn{x_1 < x_2}, where \eqn{x_1 = a + r(b-a)} and
#'    \eqn{x_2 = b - r(b-a)}, and \eqn{r=(\sqrt{5}-1)/2\approx 0.618} is the
#'    reciprocal of the golden ratio. One compares \eqn{f(x_1)} to \eqn{f(x_2)}
#'    and updates the search interval \eqn{\left[a,b\right]} as follows:
#' \itemize{
#'  \item If \eqn{f(x_1)<f(x_2)}, the solution cannot lie in
#'    \eqn{\left[x_2,b\right]};
#'    thus, update the search interval to
#'    \deqn{\left[a_\mathrm{new},b_\mathrm{new}\right]=\left[a,x_2\right]}
#'  \item If \eqn{f(x_1)>f(x_2)}, the solution cannot lie in
#'    \eqn{\left[a,x_1\right]};
#'    thus, update the search interval to
#'    \deqn{\left[a_\mathrm{new},b_\mathrm{new}\right]=\left[x_1,b\right]}
#' }
#'
#' One then chooses two new test points by replacing \eqn{a} and \eqn{b} with
#'    \eqn{a_\mathrm{new}} and \eqn{b_\mathrm{new}} and recalculating \eqn{x_1}
#'    and \eqn{x_2} based on these new endpoints. One continues iterating in
#'    this fashion until \eqn{b-a< \tau}, where \eqn{\tau} is the desired
#'    tolerance.
#'
#' @param f A function of one variable that returns a numeric value
#' @param a A numeric of length 1 representing the lower endpoint of the
#'    search interval
#' @param b A numeric of length 1 representing the upper endpoint of the
#'    search interval; must be greater than \code{a}
#' @param tol A numeric of length 1 representing the tolerance used to
#'    determine when the search interval has been narrowed sufficiently
#'    for convergence
#' @param maxitgss An integer of length 1 representing the maximum number
#'    of iterations to use in the search
#' @param ... Additional arguments to pass to \code{f}
#'
#' @return A list object containing the following:
#' \itemize{
#'  \item \code{argmin}, the argument of \code{f} that minimises \code{f}
#'  \item \code{funmin}, the minimum value of \code{f} achieved at \code{argmin}
#'  \item \code{converged}, a logical indicating whether the convergence tolerance
#'    was satisfied
#'  \item \code{iterations}, an integer indicating the number of search iterations
#'    used
#' }
#'
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[cmna]{goldsect}} is similar to this function, but does
#'    not allow the user to pass additional arguments to \code{f}
#'
#' @examples
#' f <- function(x) (x - 1) ^ 2
#' GSS(f, a = 0, b = 2)
#'

GSS <- function(f, a, b, tol = 1e-8, maxitgss = 100L, ...) {
  gtau <- (sqrt(5) - 1) / 2
  count <- 0
  if (a > b) stop("b must be greater than a")
  x <- c(a + (1 - gtau) * (b - a), a + gtau * (b - a))
  while (pracma::Norm(b - a, p = 2) > tol & count < maxitgss) {
    count <- count + 1
    f1 <- f(x[1], ...)
    f2 <- f(x[2], ...)
    if (is.na(f1) || is.nan(f1) || is.null(f1) ||
        is.na(f2) || is.nan(f2) || is.null(f2)) print(paste0("f1:", f1, "; f2: ", f2))
    if (f1 < f2) {
      b <- x[2]
      x[2] <- x[1]
      x[1] <- a + (1 - gtau) * (b - a)
    } else {
      a <- x[1]
      x[1] <- x[2]
      x[2] <- a + gtau * (b - a)
    }
  }
  if (count == maxitgss)
    warning("Golden Section Search algorithm did not converge")
  optpoint <- (a + b) / 2
  list("argmin" = optpoint, "funmin" = f(optpoint, ...),
       "converged" = pracma::Norm(b - a, p = 2) <= tol,
       "iterations" = count)
}
