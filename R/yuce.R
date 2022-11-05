#' Yüce's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the two methods of
#'    \insertCite{Yuce08;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' @details These two tests are straightforward to implement; in both cases the
#'    test statistic is a function only of the residuals of the linear
#'    regression model. The first test statistic has an asymptotic chi-squared
#'    distribution and the second has an asymptotic \eqn{t}-distribution. In
#'    both cases the degrees of freedom are \eqn{n-p}. The chi-squared test
#'    is right-tailed whereas the \eqn{t}-test is two-tailed.
#'
#' @param method A character indicating which of the two tests presented in
#'    \insertCite{Yuce08;textual}{skedastic} should be implemented. Possible
#'    values are \code{"A"} (the chi-squared test) and
#'    \code{"B"} (the \eqn{t}-test). Partial matching is used and the argument
#'    is not case-sensitive.
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
#' yuce(mtcars_lm, method = "A")
#' yuce(mtcars_lm, method = "B")
#'

yuce <- function(mainlm, method = c("A", "B"), statonly = FALSE) {

  processmainlm(m = mainlm, needy = FALSE)

  method <- match.arg(toupper(method), c("A", "B"))
  n <- length(e)

  if (method == "A") {
    alt <- "greater"
    teststat <- sum(e ^ 2) / ((sum(abs(e)) ^ 2 * pi) /
                                (pi - 2 + 2 * n * (n - p)))
    if (statonly) return(teststat)
    fullmethod <- "Chi-Squared Test"
    pval <- stats::pchisq(teststat, df = n - p, lower.tail = FALSE)
  } else if (method == "B") {
    alt <- "two.sided"
    abspart <- sum(abs(e)) ^ 2 * pi / (pi - 2 + 2 * n * (n - p))
    teststat <- (sum(e ^ 2) / (n - p) - abspart) /
      sqrt(2 / (n - p) * abspart ^ 2)
    if (statonly) return(teststat)
    fullmethod <- "t Test"
    pval <- 2 * stats::pt(abs(teststat), df = n - p, lower.tail = FALSE)
  } else stop("Invalid `method` argument")

  rval <- structure(list(statistic = teststat, p.value = pval,
               null.value = "Homoskedasticity", alternative = alt,
               method = fullmethod), class = "htest")
  broom::tidy(rval)
}
