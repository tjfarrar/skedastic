#' Harrison and McCabe's Test for Heteroskedasticity in a Linear Regression Model
#'
#' This function implements the method of
#'    \insertCite{Harrison79;textual}{skedastic} for testing for heteroskedasticity
#'    in a linear regression model.
#'
#' @details The test assumes that heteroskedasticity, if present, is
#'    monotonically related to one of the explanatory variables (known as the
#'    deflator). The OLS residuals \eqn{e} are placed in increasing order of
#'    the deflator variable and we let \eqn{A} be an \eqn{n\times n} selector
#'    matrix whose first \eqn{m} diagonal elements are 1 and all other elements
#'    are 0. The alternative hypothesis posits that the error variance changes
#'    around index \eqn{m}. Under the null hypothesis of homoskedasticity, the
#'    ratio of quadratic forms \eqn{Q=\frac{e' A e}{e'e}} should be close to
#'    \eqn{\frac{m}{n}}. Since the test statistic \eqn{Q} is a ratio of
#'    quadratic forms in the errors, the Imhof algorithm is used to compute
#'    \eqn{p}-values (with normality of errors assumed).
#'
#' @param m Either a \code{double} giving the proportion of the \eqn{n}
#'    diagonal elements of \eqn{A} that are ones, or an \code{integer}
#'    giving the index \eqn{m} up to which the diagonal elements are ones.
#'    Defaults to \code{0.5}.
#' @param alternative A character specifying the form of alternative
#'    hypothesis. If it is suspected that the error variance is positively
#'    associated with the deflator variable, \code{"less"}. If it is
#'    suspected that the error variance is negatively associated with deflator
#'    variable, \code{"greater"}. If no information is available on the suspected
#'    direction of the association, \code{"two.sided"}. Defaults to
#'    \code{"less"}.
#'
#' @inheritParams breusch_pagan
#' @inheritParams carapeto_holt
#'
#' @return An object of \code{\link[base]{class}} \code{"htest"}. If object is
#'    not assigned, its attributes are displayed in the console as a
#'    \code{\link[tibble]{tibble}} using \code{\link[broom]{tidy}}.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[lmtest:hmctest]{lmtest::hmctest}}, another
#'    implementation of the Harrison-McCabe Test. Note that the \eqn{p}-values
#'    from that function are simulated rather than computed from the
#'    distribution of a ratio of quadratic forms in normal random vectors.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' harrison_mccabe(mtcars_lm, deflator = "qsec")
#'

harrison_mccabe <- function(mainlm, deflator = NA, m = 0.5,
                    alternative = c("less", "greater", "two.sided"),
                    twosidedmethod = c("doubled", "kulinskaya"),
                    qfmethod = "imhof", statonly = FALSE) {

  alternative <- match.arg(alternative, c("less", "greater", "two.sided"))
  twosidedmethod <- match.arg(twosidedmethod, c("doubled", "kulinskaya"))

  processmainlm(m = mainlm, needy = FALSE)

  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }

  n <- nrow(X)
  checkdeflator(deflator, X, p, hasintercept[[1]])

  myA <- diag(n)
  if ((is.integer(m) || m %% 1 == 0) && (m > 1 && m < n)) {
    diag(myA)[(m + 1):n] <- 0
    m <- m / n
  } else if (is.double(m) && m %% 1 != 0 && m > 0 && m < 1) {
    mint <- round(m * n)
    if (mint %in% c(1, n)) stop("Invalid m value; when rounded, m * n must yield an integer between 2 and n - 1")
    diag(myA)[(mint + 1):n] <- 0
    m <- mint / n
  } else stop("Invalid m value; m must be an integer between 2 and n - 1 or a double between 0 and 1 exclusive")

  M <- fastM(X, n)

  if (!is.na(deflator) && !is.null(deflator)) {
    if (!is.na(suppressWarnings(as.integer(deflator)))) {
      deflator <- as.integer(deflator)
    }
    e <- e[order(X[, deflator])]
  }

  teststat <- as.double((t(e) %*% myA %*% e) / crossprod(e))
  if (statonly) return(teststat)
  if (alternative == "greater") {
    pval <- pRQF(r = teststat, A = M %*% myA %*% M,
                    B = M, algorithm = qfmethod, lower.tail = FALSE)
  } else if (alternative == "less") {
    pval <- pRQF(r = teststat, A = M %*% myA %*% M,
                    B = M, algorithm = qfmethod, lower.tail = TRUE)
  } else if (alternative == "two.sided") {
    pval <- twosidedpval(q = teststat, CDF = pRQF,
                         method = twosidedmethod, locpar = m, continuous = TRUE,
                         A = M %*% myA %*% M, B = M, algorithm = qfmethod,
                         lower.tail = TRUE)
  }

  rval <- structure(list(statistic = teststat, p.value = pval,
            parameter = m, null.value = "Homoskedasticity",
            alternative = alternative), class = "htest")
  broom::tidy(rval)
}
