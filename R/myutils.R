
columnof1s <- function (X) { # check whether a design matrix has a column of 1s
  columnsof1s <- apply(X, 2, function(x) all(x == 1))
  r <- vector("list", 2)
  r[[1]] <- any(columnsof1s)
  r[[2]] <- which(columnsof1s)
  return(r)
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

do_omit <- function(omit, n, p) {
  if (is.character(omit)) {
    omit <- match.arg(omit, c("first", "last", "random"))
    omit_passed <- omit
    if (omit == "first") {
      omit_ind <- 1:p
    } else if (omit == "last") {
      omit_ind <- (n - p + 1):n
    } else if (omit == "random") {
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

bvtnormcub <- function(x, sigma_e_i, sigma_e_j, rho) {
  function(x)
    force(1 / (2 * pi * sigma_e_i * sigma_e_j * sqrt(1 - rho ^ 2)) *
    exp(- 1 / (2 * (1 - rho ^ 2)) * ((x[1] / sigma_e_i) ^ 2 +
    (x[2] / sigma_e_j) ^ 2 - 2 * rho * x[1] * x[2]) / (sigma_e_i * sigma_e_j)))
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

supfunc <- function(B, m) {  # Function to find supremum of absolute differences in Brownian Bridge
  # as per Rackauskas and Zuokas (2007) equation (11)
  kseq <- 1:(m - 1)
  unlist(lapply(kseq, function(k)
    max(abs(B[(m - k + 1):(m + 1)] -
              B[1:(k + 1)]))))
}

rksim <- function(R, m, sqZ, seed, alpha = alpha) { # generates pseudorandom variates from T_alpha distribution
                                                    # of Rackauskas-Zuokas Test
  hseq <- (1:(m - 1)) / m
  set.seed(seed)
  T_alpha <- replicate(R, max(supfunc(B = do_brownbridge_me(stats::rnorm(m), sq = sqZ), m = m) *
                                hseq ^ (-alpha)))
}
