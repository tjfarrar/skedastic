#' Li-Yao ALRT and CVT Tests for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods of
#'    \insertCite{Li19;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' These two tests are straightforward to implement; in both cases the test
#'     statistic is a function only of the residuals of the linear regression
#'     model. Furthermore, in both cases the test statistic is asymptotically
#'     normally distributed under the null hypothesis of homoskedasticity.
#'     Both tests are right-tailed. These tests are designed to be especially
#'     powerful in high-dimensional regressions, i.e. when the number of
#'     explanatory variables is large.
#' @param method A character indicating which of the two tests derived in
#'    \insertCite{Li19;textual}{skedastic} should be implemented. Possible
#'    values are \code{"alrt"} (approximate likelihood-ratio test) and
#'    \code{"cvt"} (coefficient-of-variation test). Default is \code{"cvt"}.
#'    Partial matching is used.
#' @param baipanyin A logical. Should the central limit theorem of
#'    \insertCite{Bai16;textual}{skedastic} be used to determine the
#'    \eqn{p}-value for the coefficient-of-variation test (assuming normally
#'    distributed errors)? If \code{FALSE}, the asymptotic null distribution
#'    in \insertCite{Li19;textual}{skedastic} is used, which requires the
#'    assumption that the design variables are random and normally distributed.
#'    This argument is ignored if \code{method} is \code{"alrt"}.
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
#' li_yao(mtcars_lm, method = "alrt")
#' li_yao(mtcars_lm, method = "cvt")
#' li_yao(mtcars_lm, method = "cvt", baipanyin = FALSE)
#' # Same as first example
#' li_yao(list("e" = mtcars_lm$residuals), method = "alrt")
#'

li_yao <- function(mainlm, method = c("cvt", "alrt"), baipanyin = TRUE,
                   statonly = FALSE) {

  processmainlm(m = mainlm, needy = FALSE)

  method <- match.arg(method, c("cvt", "alrt"))
  n <- length(e)

  if (method == "alrt") {
    teststat <- log(( 1 / n * sum(e ^ 2)) / (prod(e ^ 2) ^ (1 / n)))
    if (statonly) return(teststat)
    fullmethod <- "Approximate Likelihood-Ratio Test"
    pval <- stats::pnorm((teststat - (log(2) - digamma(1))) /
                sqrt((pi ^ 2 / 2 - 2) / n), lower.tail = FALSE)
  } else if (method == "cvt") {
    mbar <- 1 / n * sum(e ^ 2)
    teststat <- (1 / n * sum((e ^ 2 - mbar) ^ 2)) / mbar ^ 2
    if (statonly) return(teststat)
    if (baipanyin) {
      M <- fastM(X, n)
      MoM <- M * M
      trMom <- sum(diag(MoM))
      a <- as.double(3 * n * trMom / ((n - p) ^ 2 + 2 * (n - p)) - 1)
      Delta <- c(n / ((n - p) ^ 2 + 2 * (n - p)), 3 * n ^ 2 * trMom /
                   ((n - p) ^ 2 + 2 * (n - p)) ^ 2)
      Theta11 <- 72 * t(diag(M)) %*% MoM %*% diag(M) + 24 * trMom ^ 2
      Theta22 <- 8 * (n - p) ^ 3 / n ^ 2
      Theta12 <- (n - p) / n * 24 * trMom
      Theta <- matrix(data = c(Theta11, Theta12, Theta12, Theta22), nrow = 2,
                      ncol = 2)
      b <- as.double(t(Delta) %*% Theta %*% Delta)
      pval <- stats::pnorm((teststat - a) / sqrt(b), lower.tail = FALSE)
    } else {
      pval <- stats::pnorm((teststat - 2) / sqrt(24 / n), lower.tail = FALSE)
    }
    fullmethod <- "Coefficient-of-Variation Test"
  } else stop("Invalid `method` argument")

  rval <- structure(list(statistic = teststat, p.value = pval,
               null.value = "Homoskedasticity", alternative = "greater",
               method = fullmethod), class = "htest")
  broom::tidy(rval)
}
