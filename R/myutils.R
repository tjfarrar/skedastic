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
      e <- stats::resid(m)
    }
    if (needy) y <- stats::model.response(stats::model.frame(m))
    if (needyhat) yhat <- stats::fitted(m)
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
      if (is.vector(X) && !is.matrix(X)) X <- as.matrix(X)
      badrows <- which(apply(cbind(y, X), 1, function(x) any(is.na(x),
                                          is.nan(x), is.infinite(x))))
      if (length(badrows) > 0) {
        message("Rows of data containing NA/NaN/Inf values removed")
        y <- y[-badrows]
        X <- X[-badrows, , drop = FALSE]
      }
    }
    if ((!exists("e", where = environment(), inherits = FALSE) || is.null(e))) {
      if (neede) {
        m <- stats::lm.fit(x = X, y = y)
        e <- stats::resid(m)
      }
      if (needyhat) yhat <- stats::fitted(m)
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
  if (is.na(d) || is.null(d)) {
  } else if (is.numeric(d) && d == as.integer(d)) {
    if (has0 && d == 1) {
      stop("deflator cannot be the model intercept")
    } else if (d > p) {
      stop("`deflator` is not the index of a column of design matrix")
    }
  } else if (is.character(d)) {
    if (d == "(Intercept)") {
      stop("deflator cannot be the model intercept")
    } else if (!(d %in% colnames(X)) &&
               suppressWarnings(is.na(as.integer(d)))) {
      stop("`deflator` is not the name of a column of design matrix")
    }
  } else stop("`deflator` must be integer or character")
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
                       ")")
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

generate_interactions <- function(X) { # generate all two-way interactions
  k <- ncol(X)
  if (k >= 2) {
    mycombn <- t(utils::combn(k, 2))
    apply(mycombn, 1, function(i, X) X[, i[1]] * X[, i[2]], X)
  } else {
    NULL
  }
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
      if (!is.na(seed.)) set.seed(seed.)
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
  X0 <- X[omit_ind, , drop = FALSE]
  X1 <- X[keep_ind, , drop = FALSE]
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

rksim <- function(R., m., sqZ., seed., alpha.) { # generates pseudorandom variates from T_alpha distribution
                                                    # of Rackauskas-Zuokas Test
  hseq <- (1:(m. - 1)) / m.
  if (!is.na(seed.)) set.seed(seed.)
  Z <- replicate(R., stats::rnorm(m.), simplify = FALSE)
  B <- lapply(X = Z, FUN = do_brownbridge_me, sq = sqZ.)
  vapply(X = B, FUN = function(b) max(supfunc(B = b, m.. = m.) *
                      hseq ^ (-alpha.)), NA_real_)
}

gqind <- function(n, p, propremove, propgroup1) {
 nremove <- as.integer(round(n * propremove))
 k <- nremove / n
 nprime <- n - nremove
 n1prime <- as.integer(round(nprime * propgroup1))
 w <- n1prime / nprime
 n2prime <- nprime - n1prime

 if (min(n1prime, n2prime) <= p) {
   maxw <- ifelse(n %% 2 == 0, 0.5, (n - 1) / (2 * n))
   n1prime <- n * maxw
   n2prime <- n - n1prime

   if (min(n1prime, n2prime) <= p) {
     stop("Your choice of `prop_central` and `group1prop` results in a subset containing less than or equal to p observations, where p is the number of coefficients to be estimated. As a result, the models cannot be fit. Even if `prop_central` were set to 0 and `propgroup1` to 0.5, the number of observations in at least one of the groups would be <= p. Thus, the Goldfeld-Quandt Test cannot be implemented.")
   } else {
     message(paste0("Your choice of `prop_central` and `group1prop` results in a subset containing less than or equal to p observations, where p is the number of coefficients to be estimated. As a result, the models cannot be fit. The test has instead been implemented with prop_central set to 0 and group1prop to ", maxw))
   }
 }
 return(list(1:n1prime, (n - n2prime + 1):n))
}

fastM <- function(X, n) {
  diag(n) - X %*% solve(crossprod(X)) %*% t(X)
}

SKH <- function(i) {
  2 * (1 - cos((pi * i) / (length(i) + 1)))
}

is.btwn01 <- function(x) {
  if (is.finite(x)) {
    (x >= 0 && x <= 1)
  } else {
    FALSE
  }
}

is.square.mat <- function(a) {
  if (!is.matrix(a))
    stop("argument a is not a matrix")
  return(nrow(a) == ncol(a))
}

is.symmetric.mat <- function(a) {
  if (!is.matrix(a)) {
    stop("argument a is not a matrix")
  }
  if (!is.numeric(a)) {
    stop("argument a is not a numeric matrix")
  }
  if (!is.square.mat(a))
    stop("argument a is not a square numeric matrix")

  all(berryFunctions::almost.equal(a[upper.tri(a)], t(a)[upper.tri(t(a))]))
}

is.pos.semidef.mat <- function(a) {
  if (!is.square.mat(a))
    stop("argument a is not a square matrix")
  if (!is.symmetric.mat(a))
    stop("argument a is not a symmetric matrix")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix")
  myeigen <- eigen(a, only.values = TRUE)$values
  cplxeigen <- which(Im(myeigen) != 0)

  if (length(cplxeigen) == 0 || all(berryFunctions::almost.equal(Im(myeigen[cplxeigen]), 0))) {
    myeigen[cplxeigen] <- Re(myeigen[cplxeigen])
    myeigen <- as.numeric(myeigen)
  } else stop("matrix a has complex eigenvalues")

  negeigen <- which(myeigen < 0)
  (length(negeigen) == 0 | all(berryFunctions::almost.equal(myeigen[negeigen],
                                                            0)))
}

is.singular.mat <- function(a) {
  if (!is.square.mat(a))
    stop("argument a is not a square matrix")
  if (!is.numeric(a))
    stop("argument a is not a numeric matrix")
  berryFunctions::almost.equal(det(a), 0)
}

errormat <- function(res, design, HCCME = "HC4", k = 0.7, gamma = c(1, 1.5),
                     coeff = FALSE) {
  n <- nrow(design)
  p <- ncol(design)
  if (HCCME == "const") {
    sum(res ^ 2) / (n - p) * diag(n)
  } else if (HCCME == "HC0") {
    diag(res ^ 2)
  } else if (HCCME == "HC1") {
    diag(res ^ 2 * n / (n - p))
  } else if (HCCME %in% c("HC2", "HC3", "HC4", "HC5", "HC4m")) {
    hdiag <- diag(design %*% solve(crossprod(design)) %*% t(design))
    hbar <- p / n
    if (HCCME == "HC2") {
      delta <- 1
    } else if (HCCME == "HC3") {
      delta <- 2
    } else if (HCCME == "HC4") {
      delta <- vapply(1:n, function(i) {
        min(4, hdiag[i] / hbar)
      }, NA_real_)
    } else if (HCCME == "HC5") {
      hmax <- max(hdiag)
      maxcompare <- max(4, k * hmax / hbar)
      delta <- vapply(1:n, function(i) {
        min(maxcompare, hdiag[i] / hbar)
      }, NA_real_) / 2
    } else if (HCCME == "HC4m") {
      delta <- vapply(1:n, function(i) {
        min(gamma[1], hdiag[i] / hbar) + min(gamma[2], hdiag[i] / hbar)
      }, NA_real_)
    }
    if (!coeff) {
      diag(res ^ 2 / (1 - hdiag) ^ delta)
    } else {
      solve(crossprod(design)) %*% t(design) %*%
        diag(res ^ 2 / (1 - hdiag) ^ delta) %*%
        design %*% solve(crossprod(design))
    }
  }
}

wma <- function(x, order.by2 = NULL, k2, wt2, ntrim2) {
  n <- length(x)
  if (!is.null(order.by2)) {
    o <- order(order.by2)
    x <- x[o]
  }
  w <- c(1:(k2 + 1), k2:1)
  avg <- vapply((k2 + 1):(n - k2), function(i) {
    if (wt2 == "const") {
      mean(x[(i - k2):(i + k2)], trim = ntrim2 / (2 * k2 + 1))
    } else if (wt2 == "integer") {
      sum(w * x[(i - k2):(i + k2)]) / (k2 + 1) ^ 2
    }
  }, NA_real_)
  xnew <- c(rep(NA_real_, k2), avg, rep(NA_real_, k2))
  xnew[order(o)]
}

midpoint <- function(x) {
  d <- diff(x)
  xlesslast <- x[-length(x)]
  xlesslast + d / 2
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
