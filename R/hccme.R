#' Heteroskedasticity-Consistent Covariance Matrix Estimators for
#'    Linear Regression Models
#'
#' Computes an estimate of the \eqn{n\times n} covariance matrix \eqn{\Omega}
#'    (assumed to be diagonal) of the random error vector of a linear
#'    regression model, using a specified method
#'
#' @param hcnum A character corresponding to a subscript in the name of an
#'    HCCME according to the usual nomenclature \eqn{\mathrm{HC\#}}.
#'    Possible values are:
#'    \itemize{
#'     \item "3", the default, corresponding to HC3
#'    \insertCite{MacKinnon85}{skedastic}
#'     \item "0", corresponding to HC0 \insertCite{White80}{skedastic}
#'     \item "1", corresponding to HC1 \insertCite{MacKinnon85}{skedastic}
#'     \item "2", corresponding to HC1 \insertCite{MacKinnon85}{skedastic}
#'     \item "4", corresponding to HC4 \insertCite{Cribari04}{skedastic}
#'     \item "5", corresponding to HC5 \insertCite{Cribari07}{skedastic}
#'     \item "6", corresponding to HC6 \insertCite{Aftab16}{skedastic}
#'     \item "7", corresponding to HC7 \insertCite{Aftab18}{skedastic}
#'     \item "4m", corresponding to HC4m \insertCite{Cribari11}{skedastic}
#'     \item "5m", corresponding to HC5m \insertCite{Li17}{skedastic}
#'     \item "const", corresponding to the homoskedastic estimator,
#'     \eqn{(n-p)^{-1}\displaystyle\sum_{i=1}^{n}e_i^2}
#'    }
#' @param sandwich A logical, defaulting to \code{FALSE}, indicating
#'    whether or not the sandwich estimator
#'    \deqn{\mathrm{Cov}{\hat{\beta}}=(X'X)^{-1}X'\hat{\Omega}X(X'X)^{-1}}
#'    should be returned instead of \eqn{\mathrm{Cov}(\epsilon)=\hat{\Omega}}
#' @param as_matrix A logical, defaulting to \code{TRUE}, indicating whether
#'    a covariance matrix estimate should be returned rather
#'    than a vector of variance estimates
#' @inheritParams breusch_pagan
#'
#' @return A numeric matrix (if \code{as_matrix} is \code{TRUE}) or else a
#'    numeric vector
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[sandwich]{vcovHC}}
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' Omega_hat <- hccme(mtcars_lm, hcnum = "4")
#' Cov_beta_hat <- hccme(mtcars_lm, hcnum = "4", sandwich = TRUE)
#'

hccme <- function(mainlm,
                  hcnum = c("3", "0", "1", "2", "4", "5",
                            "6", "7", "4m", "5m", "const"),
                  sandwich = FALSE,
                  as_matrix = TRUE) {

  processmainlm(m = mainlm)
  esq <- e ^ 2
  n <- length(esq)

  hcnum <- match.arg(hcnum,
                     c("3", "0", "1", "2", "4", "5",
                       "6", "7", "4m", "5m", "const"))

  H <- X %*% Rfast::spdinv(crossprod(X)) %*% t(X)
  h <- diag(H)
  if (hcnum %in% c("4", "5", "4m", "5m", "7")) {
    hbar <- mean(h)
    hmax <- max(h)
  }
  Omegahat <- if (hcnum == "0") {
    esq
  } else if (hcnum == "1") {
    esq * n / (n - p)
  } else if (hcnum == "2") {
    esq / (1 - h) ^ 1
  } else if (hcnum == "3") {
    esq / (1 - h) ^ 2
  } else if (hcnum == "4") {
    esq / (1 - h) ^ (min(4, h / hbar))
  } else if (hcnum == "5") {
    esq / (1 - h) ^ (1/2 * min(h / hbar, max(4, 0.7 * hmax)))
  } else if (hcnum == "4m") {
    esq / (1 - h) ^ (min(1, h / hbar) + min(1.5, h / hbar))
  } else if (hcnum == "5m") {
    esq / (1 - h) ^ (min(1, h / hbar) + min(h / hbar, max(4, 0.7 * hmax / hbar)))
  } else if (hcnum == "6") {
    omegahat <- sum(esq) / (n - p)
    cooksd <- esq * h / (omegahat * p * (1 - h) ^ 2)
    esq * sqrt(cooksd)
  } else if (hcnum == "7") {
    esq / (1 - h) ^ (min(h / hbar, sqrt(hmax / (2 * hbar))))
  } else if (hcnum == "const") {
    rep(sum(esq) / (n - p), n)
  }

  if (sandwich) {
    XtXinv <- Rfast::spdinv(crossprod(X))
    result <- XtXinv %*% t(X) %*% diag(Omegahat) %*% X %*% XtXinv
    if (as_matrix) {
      result
    } else {
      unname(diag(result))
    }
  } else {
    if (as_matrix) {
      diag(Omegahat)
    } else {
      unname(Omegahat)
    }
  }
}
