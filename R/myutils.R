processmainlm <- function(m, needX = TRUE, needy = TRUE, neede = TRUE,
                          needyhat = FALSE, needp = TRUE) {
# process mainlm object passed to function

  y <- NULL
  X <- NULL
  yhat <- NULL
  p <- NULL
  if (class(m) == "lm") {
    if (needX) X <- stats::model.matrix(m)
    if (!exists("e", where = environment(), inherits = FALSE) && neede) {
      e <- m$residuals
    }
    if (needy) y <- stats::model.response(stats::model.frame(m))
    if (needyhat) yhat <- m$fitted.values
  } else if (class(m) == "list") {
    if (is.null(names(m))) {
      y <- m[[1]]
      X <- m[[2]]
      if (length(m) > 2) e <- m[[3]]
    } else {
      y <- m$y
      X <- m$X
      e <- m$e
    }
    if (!is.null(X)) {
      badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x),
                                          is.nan(x), is.infinite(x))))
      if (length(badrows) > 0) {
        warning("Rows of data containing NA/NaN/Inf values removed")
        y <- y[-badrows]
        X <- X[-badrows, drop = FALSE]
      }
    }
    if ((!exists("e", where = environment(), inherits = FALSE) || is.null(e))) {
      if (neede) {
        m <- stats::lm.fit(X, y)
        e <- m$residuals
      }
      if (needyhat) yhat <- m$fitted.values
    } else {
      if (needyhat) yhat <- y - e
    }
  }
  if (needp) p <- ncol(X)
  if (exists("e", where = environment(), inherits = FALSE)) {
    invisible(list2env(x = list("y" = y, "X" = X, "e" = e,
                     "yhat" = yhat, "p" = p), envir = parent.frame()))
  } else {
    invisible(list2env(x = list("y" = y, "X" = X, "yhat" = yhat,
                    "p" = p), envir = parent.frame()))
  }
}

columnof1s <- function (X) { # check whether a design matrix has a column of 1s
  columnsof1s <- apply(X, 2, function(x) all(x == 1))
  r <- vector("list", 2)
  r[[1]] <- any(columnsof1s)
  r[[2]] <- which(columnsof1s)
  return(r)
}

checkdeflator <- function(d, X, p, has0) {
  if (is.numeric(d) && d == as.integer(d)) {
    if (has0 && d == 1) {
      stop("deflator cannot be the model intercept")
    } else if (d > p) {
      stop("`deflator` is not the index of a column of design matrix")
    }
  } else if (is.character(d)) {
    if (d == "(Intercept)") {
      stop("deflator cannot be the model intercept")
    } else if (!d %in% colnames(X)) {
      stop("`deflator` is not the name of a column of design matrix")
    }
  } else if (!is.null(d)) stop("`deflator` must be integer or character")
}

plotdim <- function(x) { # find all positive factors of a positive integer
  x <- as.integer(x)
  if (x < 0) stop("Positive integers only")
  div <- seq_len(x)
  factors <- div[x %% div == 0L]
  twofactors <- cbind(factors, x / factors, deparse.level = 0)
  twofactors <- twofactors[twofactors[, 1] <= twofactors[, 2], , drop = FALSE]
  twofactors[rowSums(twofactors) == min(rowSums(twofactors)), ]
}

parselabels <- function(labname, theaxis) {

  if (theaxis == "y") {
    v <- unlist(lapply(labname, function(char) substr(char, start = 1,
                      stop = max(gregexpr("_", char)[[1]]) - 1)))
    f <- unlist(lapply(labname, function(char) substr(char,
                      start = max(gregexpr("_", char)[[1]]) + 1,
                      stop = nchar(char))))
    innerexpr <- switch(v, "res" = "e[i]",
                        "res_blus" = "tilde(e)[i]", "res_stand" = "frac(e[i], hat(sigma))",
                        "res_constvar" = "frac(e[i], sqrt(m[ii]))",
                        "res_stud" = "frac(e[i], sqrt(s ^ 2 * m[ii]))",
                        v)
    if (!suppressWarnings(is.na(as.integer(f)))) {
      prefix <- switch(v, "res" = "",
                       "res_blus" = "", "res_stand" = "bgroup(\"(\",",
                       "res_constvar" = "bgroup(\"(\",",
                       "res_stud" = "bgroup(\"(\",",
                       v)
      suffix <- switch(v, "res" = "",
                       "res_blus" = "", "res_stand" = ",\")\")",
                       "res_constvar" = ",\")\")",
                       "res_stud" = ",\")\")",
                       v)
      parse(text = paste0(prefix, innerexpr, suffix, "^", f))
    } else {
      prefix <- switch(f, "sqrt" = "sqrt(", "identity" = "", "log" = "log(",
                       "abs" = "abs(phantom(i)*",
                       paste0(f, "("))
      suffix <- switch(f, "sqrt" = ")", "identity" = "", "log" = ")",
                       "abs" = "*phantom(i))",
                       paste0(f, "("))
      parse(text = paste0(prefix, innerexpr, suffix))
    }
  } else if (theaxis == "x") {
    texttoparse <- switch(labname, "index" = "i", "fitted.values" = "hat(y)[i]",
                          "fitted.values2" = "m[ii] * hat(y)[i]", labname)
    if (substr(texttoparse, 1, 1) == "X" &&
        suppressWarnings(!is.na(as.integer(substr(texttoparse, 2, nchar(texttoparse)))))) {
         j <- as.integer(substr(texttoparse, 2, nchar(texttoparse)))
         texttoparse <- paste0("X[i * ", j, "]")
    } else if (substr(texttoparse, 1, 4) == "logX" &&
      suppressWarnings(!is.na(as.integer(substr(texttoparse, 5, nchar(texttoparse)))))) {
        j <- as.integer(substr(texttoparse, 5, nchar(texttoparse)))
        texttoparse <- paste0("log(X[i * ", j, "])")
    }
    parse(text = texttoparse)
  }
}

generate_interactions <- function (X) { # generate all two-way interactions
  k <- ncol(X)
  mycombn <- t(utils::combn(k, 2))
  apply(mycombn, 1, function(i, X) X[, i[1]] * X[, i[2]], X)
}

checklm <- function(mylm) { # Check validity of lm object
  if (class(mylm) != "lm") stop("`mylm` must be an object of class `lm`")
  if (!all(c("residuals", "df.residual", "coefficients", "model",
             "fitted.values", "terms") %in% names(mylm))) {
    stop("`mylm` must contain the components expected in class `lm`")
  }
  if (any(is.na(stats::model.frame(mylm)))) {
    stop("Missing values in `model.frame` not allowed.")
  }
}

margin_indices <- function(v, p, sub_ind) {

  k <- length(v)
  reptimes <- floor(p / (k - 1))
  pmodk1 <- p %% (k - 1)
  omit_ind <- vector("list", reptimes + 1)
  if (reptimes >= 1) {
    for (r in 1:reptimes) {
      if (r %% 2 == 1) {
        omit_ind[[r]] <- unlist(lapply(sub_ind, max))[1:(k - 1)]
      } else if (r %% 2 == 0) {
        omit_ind[[r]] <- unlist(lapply(sub_ind, min))[2:k]
      }
    }
  }
  if (pmodk1 != 0) {
    if (reptimes %% 2 == 0) {
      omit_ind[[reptimes + 1]] <- unlist(lapply(sub_ind, max))[1:pmodk1]
    } else if (reptimes %% 2 == 1) {
      omit_ind[[reptimes + 1]] <- unlist(lapply(sub_ind, min))[2:(pmodk1 + 1)]
    }
  }
  unlist(omit_ind)
}

do_omit <- function(omit, n, p, seed.) {
  if (is.character(omit)) {
    omit <- match.arg(omit, c("first", "last", "random"))
    omit_passed <- omit
    if (omit == "first") {
      omit_ind <- 1:p
    } else if (omit == "last") {
      omit_ind <- (n - p + 1):n
    } else if (omit == "random") {
      if (!is.null(seed.)) set.seed(seed.)
      omit_ind <- sort(sample(1:n, p))
    } else stop("Invalid character for `omit` argument")
  } else if (is.numeric(omit)) {
    if (length(omit) != p) stop("length of `omit` does not equal number of
                                columns of original design matrix")
    omit_passed <- "intvector"
    omit_ind <- as.integer(omit)
    if (is.unsorted(omit_ind)) omit_ind <- sort(omit_ind)
    if (min(omit_ind) < 1 || max(omit_ind) > n) stop("`omit` contains
                                                     subscripts out of range")
  } else stop("Invalid `omit` argument: must be numeric vector or character
              value")

  return(list("omit_passed" = omit_passed, "omit_ind" = omit_ind))
}

do_Xmats <- function(X, n, p, omit_ind) {
  keep_ind <- setdiff(1:n, omit_ind)
  X0 <- X[omit_ind, ]
  X1 <- X[keep_ind, ]
  if (identical(omit_ind, 1:p)) {
    X_ord <- X
  } else {
    X_ord <- rbind(X0, X1)
  }
  X_ord_sq <- t(X_ord) %*% X_ord
  return(list("X0" = X0, "X1" = X1, "X_ord_sq" = X_ord_sq))
}

normexpect_integrand <- function(x, sigma1, sigma2, rho) {
  function(x)
    force(sqrt(abs(x[1] * x[2])) / (2 * pi * sigma1 * sigma2 * sqrt(1 - rho ^ 2)) *
    exp(- 1 / (2 * (1 - rho ^ 2)) * ((x[1] / sigma1) ^ 2 +
    (x[2] / sigma2) ^ 2 - 2 * rho * x[1] * x[2] / (sigma1 * sigma2))))
}

value_possible <- function(x, myCDF, ...) {
  if (x %% 1 == 0) {
    if (myCDF(x, ...) - myCDF(x - 1, ...) == 0) {
      FALSE
    } else if (myCDF(x, ...) - myCDF(x - 1, ...) > 0) {
      TRUE
    }
  } else if (x %% 1 > 0) {
    FALSE
  }
}

meanfromCDF <- function(theCDF, cont, suplim, ...) {
  CDF2 <- function(x) theCDF(x, ...)
  surv <- function(x) 1 - theCDF(x, ...)
  if (cont) {
    stats::integrate(surv, lower = 0, upper = Inf)$value -
      stats::integrate(CDF2, lower = -Inf, upper = 0)$value
  } else {
    if (missing(suplim)) {
      meanval <- sum(surv(0:1e6)) - sum(CDF2(-1e6:-1))
    } else {
      if (is.infinite(suplim[2])) {
        suplim[2] <- 1e6
        warning("Infinite upper limit of support truncated at 1e6")
      }
      if (is.infinite(suplim[1])) {
        suplim[1] <- -1e6
        warning("Infinite lower limit of support truncated at -1e6")
      }
      negsupport <- suplim[1] < 0
      possupport <- suplim[2] > 0
      if (negsupport && possupport) {
        meanval <- sum(surv(0:suplim[2])) - sum(CDF2(suplim[1]:-1))
      } else if (negsupport && !possupport) {
        meanval <- 0 - sum(CDF2(suplim[1]:suplim[2]))
      } else if (!negsupport && possupport) {
        meanval <- sum(surv(suplim[1]:suplim[2]))
      }
    }
    if (meanval %% 1 < 1e-4 || meanval %% 1 > (1 - 1e-4)) {
      as.integer(round(meanval))
    } else {
      meanval
    }
  }
}

# normexp_integrand2 <- function(x, Sigmat) {
#   function(x)
#     force(mvtnorm::dmvnorm(x, sigma = Sigmat))
# }

# qreg arguments:
# x and y are independent and dependent variables (matrix and vector)
# qval is the quantile to be used, which defaults to median
# the argument q can be used instead of qval to alter the quantile used
# the argument pr ???
# the argument xout indicates whether data are to be removed if
#   the value of x is flagged as an outlier
# the argument outfun indicates the method to detect outliers, which defaults
#   to the MAD-median rule when there is only one independent variable
# plotit determines whether to draw a scatter plot
# xlab and ylab would be axis labels of plot
# op = 1, v2 = TRUE, method = "br", WARN = FALSE ???


do_brownbridge_me <- function(z, sq = FALSE) { # Function to transform independent std normal variates into Brownian Bridge
  m <- length(z)
  if (sq) {
    wiener <- stats::ts(cumsum(z ^ 2 / sqrt(2 * m)), start = 1 / m, frequency = m)
    stats::ts(c(0, wiener - stats::time(wiener) * as.vector(wiener)[m]),
        start = 0, frequency = m)
  } else {
    wiener <- stats::ts(cumsum(z / sqrt(m)), start = 1 / m, frequency = m)
    stats::ts(c(0, wiener - stats::time(wiener) * as.vector(wiener)[m]),
       start = 0, frequency = m)
  }
}

supfunc <- function(B, m..) {  # Function to find supremum of absolute differences in Brownian Bridge
  # as per Rackauskas and Zuokas (2007) equation (11)
  kseq <- 1:(m.. - 1)
  unlist(lapply(kseq, function(k)
    max(abs(B[(m.. - k + 1):(m.. + 1)] -
              B[1:(k + 1)]))))
}

rksim <- function(R., m., sqZ., seed., alpha. = alpha) { # generates pseudorandom variates from T_alpha distribution
                                                    # of Rackauskas-Zuokas Test
  hseq <- (1:(m. - 1)) / m.
  if (!is.null(seed.)) set.seed(seed.)
  Z <- replicate(R., stats::rnorm(m.), simplify = FALSE)
  B <- lapply(X = Z, FUN = do_brownbridge_me, sq = sqZ.)
  vapply(X = B, FUN = function(b) max(supfunc(B = b, m.. = m.) *
                      hseq ^ (-alpha.)), NA_real_)
}

gqind <- function(n, p, g) {
 k <- as.integer(round(n * p))
 nprime <- n - k
 n1prime <- as.integer(round(nprime * g))
 n2prime <- nprime - n1prime
 if (min(n1prime, n2prime) <= p) {
   stop("Your choice of `prop_central` and `group1prop` results in a subset
        containing less than or equal to p observations, where p is the number
        of coefficients to be estimated. As a result, the models cannot be fit.")
 } else if (k < 0) {
   stop("`prop_central` cannot be negative")
 }
 return(list(1:n1prime, (n - n2prime + 1):n))
}

fastM <- function(X, n) {
  diag(n) - X %*% solve(crossprod(X)) %*% t(X)
}

IRfunc <- function(X, e, B, method) {
  p <- ncol(X)
  n <- length(e)
  sigma_hat_sq <- sum(e ^ 2) / n
  H <- X %*% solve(crossprod(X)) %*% t(X)
  xi <- replicate(B, stats::rnorm(n))
  if (method == "Xj") {
    Hminus <- lapply(1:p, function(j) X[, -j, drop = FALSE] %*%
                       solve(t(X[, -j, drop = FALSE]) %*% X[, -j, drop = FALSE]) %*%
                       t(X[, -j, drop = FALSE]))
    w <- sapply(1:p, function(j) diag(H) - diag(Hminus[[j]]))
    IR <- vapply(1:p, function(j) sum(w[, j] * e ^ 2) / sigma_hat_sq, NA_real_)
    W <- sqrt(n) * (IR - 1)
    Wstar <- sapply(1:B, function(b) vapply(1:p, function(j) sqrt(n) *
            sum(xi[, b] * ((w[, j] - 1 / n * IR[j]) * (e ^ 2 / sigma_hat_sq - 1) -
            1 / n * (IR[j] - 1))), NA_real_))
    Px <- vapply(1:3, function(j) sum(Wstar[j, ] >= W[j]) / B, NA_real_)
    return(list(W, Px))
  } else if (method == "pool") {
    wpool <- diag(H) / p
    IRpool <- sum(e ^ 2 * wpool) / sigma_hat_sq
    Wstarpool <- vapply(1:B, function(b) sqrt(n) *
                          sum(xi[, b] * ((wpool - 1 / n * IRpool) *
                                           (e ^ 2 / sigma_hat_sq - 1) - 1 / n * (IRpool - 1))), NA_real_)
    Wpool <- sqrt(n) * (IRpool - 1)
    P <- sum(Wstarpool >= Wpool) / B
    return(list(Wpool, P))
  }
}

SKH <- function(i) {
  2 * (1 - cos((pi * i) / (length(i) + 1)))
}
