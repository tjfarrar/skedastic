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
#'    values are "alrt" (approximate likelihood-ratio test) and "cvt"
#'    (coefficient-of-variation test). Default is "alrt". It is acceptable to
#'    specify only the first letter.
#' @inheritParams breusch_pagan
#'
#' @return An object of \code{\link[base]{class}} "htest". If object is not
#'    assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' li_yao(mtcars_lm, method = "alrt")
#' li_yao(mtcars_lm, method = "cvt")
#'

li_yao <- function (mainlm, method = "alrt") {

  if (class(mainlm) == "list") {
    y <- mainlm[[1]]
    X <- mainlm[[2]]
    badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x), is.nan(x), is.infinite(x))))
    if (length(badrows) > 0) {
      warning("Rows of data containing NA/NaN/Inf values removed")
      y <- y[-badrows]
      X <- X[-badrows, ]
    }
    mainlm <- stats::lm.fit(X, y)
  }

  method <- match.arg(method, c("alrt", "cvt"))
  e <- mainlm$residuals
  n <- length(e)

  if (method == "alrt") {
    T1 <- log(( 1 / n * sum(e ^ 2)) / (prod(e ^ 2) ^ (1 / n)))
    teststat <- (T1 - (log(2) - digamma(1))) / sqrt((pi ^ 2 / 2 - 2) / n)
    fullmethod <- "Approximate Likelihood-Ratio Test"
  } else if (method == "cvt") {
    mbar <- 1 / n * sum(e ^ 2)
    T2 <- (1 / n * sum((e ^ 2 - mbar) ^ 2)) / mbar ^ 2
    teststat <- (T2 - 2) / sqrt(24 / n)
    fullmethod <- "Coefficient-of-Variation Test"
  } else stop("Invalid `method` argument")

  pval <- 1 - stats::pnorm(teststat)
  rval <- structure(list(statistic = teststat, p.value = pval,
               null.value = "Homoskedasticity", alternative = "Heteroskedasticity",
               method = fullmethod), class = "htest")
  broom::tidy(rval)
}
