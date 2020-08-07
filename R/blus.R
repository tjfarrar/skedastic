#' Compute Best Linear Unbiased Scalar-Covariance (BLUS) residuals from a linear model
#'
#' This function computes the Best Linear Unbiased Scalar-Covariance (BLUS)
#'    residuals from a linear model, as defined in
#'    \insertCite{Theil65;textual}{skedastic} and explained further in
#'    \insertCite{Theil68;textual}{skedastic}.
#'
#'    Under the ideal linear model conditions, the BLUS residuals have a scalar
#'    covariance matrix \eqn{\sigma^2 I} (meaning they have a constant variance
#'    and are mutually uncorrelated), unlike the OLS residuals, which have
#'    covariance matrix \eqn{\sigma^2 M} where \eqn{M} is a function of the
#'    design matrix. Use of BLUS residuals could improve the performance of
#'    tests for heteroskedasticity and/or autocorrelation in the linear model.
#'    A linear model with \eqn{n} observations and an \eqn{n\times p} design
#'    matrix yields only \eqn{n-p} BLUS residuals. The choice of which \eqn{p}
#'    observations will not be represented in the BLUS residuals is specified
#'    within the algorithm.
#'
#' @param omit A numeric vector of length \eqn{p} (the number of columns in the
#'    linear model design matrix) giving the indices of \eqn{p} observations to omit in
#'    the BLUS residual vector; or a character partially matching \code{"first"}
#'    (for the first \eqn{p}) observations, \code{"last"} (for the last \eqn{p}
#'    observations), or \code{"random"} (for a random sample of \eqn{p} indices
#'    between 1 and \eqn{n}). Defaults to \code{"first"}.
#' @param keepNA A logical. Should BLUS residuals for omitted observations be
#'    returned as \code{NA_real_} to preserve the length of the residual vector?
#' @param exhaust An integer. If singular matrices are encountered
#'    using the passed value of \code{omit}, how many random combinations
#'    of \eqn{p} indices should be attempted before an error is thrown? If
#'    \code{NA} (the default), all possible combinations are attempted
#'    provided that \eqn{{n \choose p} \le 10^4}; otherwise up to
#'    \eqn{10^4} random samples of size \code{p} from \code{1:n} are
#'    attempted (with replacement). Integer values of \code{exhaust}
#'    greater than \code{1e4L} are treated as \code{NA}.
#' @inheritParams breusch_pagan
#' @inheritParams wilcox_keselman
#'
#' @return A double vector of length \eqn{n} containing the BLUS residuals
#'    (with \code{NA_real_}) for omitted observations), or a double vector
#'    of length \eqn{n-p} containing the BLUS residuals only (if \code{keepNA}
#'    is set to \code{FALSE})
#'
#' @references{\insertAllCited{}}
#' @export
#' @seealso H. D. Vinod's online article,
#'    \href{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2412740}{Theil's
#'    BLUS Residuals and R Tools for Testing and Removing Autocorrelation and
#'    Heteroscedasticity}, for an alternative function for computing BLUS
#'    residuals.
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' blus(mtcars_lm)
#' plot(mtcars_lm$residuals, blus(mtcars_lm))
#' # Same as first example
#' mtcars_list <- list("y" = mtcars$mpg, "X" = cbind(1, mtcars$wt, mtcars$qsec, mtcars$am))
#' blus(mtcars_list)
#' # Again same as first example
#' mtcars_list2 <- list("e" = mtcars_lm$residuals, "X" = cbind(1, mtcars$wt, mtcars$qsec, mtcars$am))
#' blus(mtcars_list2)
#' # BLUS residuals cannot be computed with `omit = "last"` in this example, so
#' # omitted indices are randomised:
#' blus(mtcars_lm, omit = "last")
#'

blus <- function(mainlm, omit = c("first", "last", "random"), keepNA = TRUE,
                  exhaust = NA, seed = 1234) {

  processmainlm(m = mainlm, needy = FALSE)

  n <- nrow(X)
  if (!is.na(seed)) set.seed(seed)
  omitfunc <- do_omit(omit, n, p, seed)
  Xmats <- do_Xmats(X, n, p, omitfunc$omit_ind)
  singular_matrix <- FALSE

  if (is.singular.mat(Xmats$X_ord_sq) ||
      is.singular.mat(Xmats$X0)) {

    singular_matrix <- TRUE
    message("Passed `omit` argument resulted in singular matrix; BLUS residuals
          could not be computed. Randomly chosen combinations of indices to
          omit will be attempted according to `exhaust` argument passed.")
    ncombn <- choose(n, p)
    if ((is.na(exhaust) || is.null(exhaust))) {
      dosample <- (ncombn > 1e4)
      numsample <- 1e4
    } else if (exhaust <= 0) {
      stop("`exhaust` is not positive; no attempts will be made to find subset to omit")
    } else {
      dosample <- (ncombn > exhaust)
      numsample <- min(ncombn, exhaust)
    }

    if (dosample) {
      subsetstotry <- unique(t(replicate(numsample, sort(sample(x = n, size = p)))))
      rowstodo <- 1:nrow(subsetstotry)
    } else {
      subsetstotry <- t(utils::combn(n, p))
      maxrow <- min(exhaust, nrow(subsetstotry), na.rm = TRUE)
      rowstodo <- sample(1:nrow(subsetstotry), maxrow, replace = FALSE)
    }

    for (r in rowstodo) {
      omitfunc <- do_omit(subsetstotry[rowstodo[r], , drop = FALSE], n, p)
      Xmats <- do_Xmats(X, n, p, omitfunc$omit_ind)
      if (!is.singular.mat(Xmats$X_ord_sq) &&
          !is.singular.mat(Xmats$X0)) {
        singular_matrix <- FALSE
        message(paste0("Success! Subset of indices found that does not yield singular
                   matrix: ", paste(omitfunc$omit_ind, collapse = ",")))
        break
      }
    }
  }
  if (singular_matrix) stop("No subset of indices to omit was found that
                            avoided a singular matrix.")

  keep_ind <- setdiff(1:n, omitfunc$omit_ind)
  G <- Xmats$X0 %*% solve(Xmats$X_ord_sq) %*% t(Xmats$X0)
  Geig <- eigen(G, symmetric = TRUE)
  lambda <- sqrt(Geig$values)
  q <- as.data.frame(Geig$vectors)
  Z <- Reduce(`+`, mapply(`*`, lambda / (1 + lambda),
                          lapply(q, function(x) tcrossprod(x)),
                          SIMPLIFY = FALSE))
  e0 <- e[omitfunc$omit_ind]
  e1 <- e[keep_ind]
  e_tilde <- c(e1 - Xmats$X1 %*% solve(Xmats$X0) %*% Z %*% e0)

  if (keepNA) {
    rval <- rep(NA_real_, n)
    rval[keep_ind] <- e_tilde
  } else {
    rval <- e_tilde
  }
  rval
}
