#' Auxiliary Nonlinear Variance Model
#'
#' Fits an Auxiliary Nonlinear Variance Model (ANLVM) to estimate the error
#'    variances of a heteroskedastic linear regression model.
#'
#' @details The ANLVM model equation is
#'    \deqn{e_i^2=\displaystyle\sum_{k=1}^{n} g(X_{k\cdot}'\gamma) m_{ik}^2+u_i},
#'    where \eqn{e_i} is the \eqn{i}th Ordinary Least Squares residual,
#'    \eqn{X_{k\cdot}} is a vector corresponding to the \eqn{k}th row of the
#'    \eqn{n\times p} design matrix \eqn{X}, \eqn{m_{ik}^2} is the
#'    \eqn{(i,k)}th element of the annihilator matrix \eqn{M=I-X(X'X)^{-1}X'},
#'    \eqn{u_i} is a random error term, \eqn{\gamma} is a \eqn{p}-vector of
#'    unknown parameters, and \eqn{g(\cdot)} is a continuous, differentiable
#'    function that need not be linear in \eqn{\gamma}, but must be expressible
#'    as a function of the linear predictor \eqn{X_{k\cdot}'\gamma}.
#'    This method has been developed as part of the author's doctoral research
#'    project.
#'
#' The parameter vector \eqn{\gamma} is estimated using the maximum
#'    quasi-likelihood method as described in section 2.3 of
#'    \insertCite{Seber03;textual}{skedastic}. The optimisation problem is
#'    solved numerically using a Gauss-Newton algorithm.
#'
#' For further discussion of feature selection and the methods for choosing the
#'    number of clusters to use with the clustering version of the model, see
#'    \code{\link{alvm.fit}}.
#'
#' @param g A numeric-valued function of one variable, or a character denoting
#'    the name of such a function. \code{"sq"} is allowed as a way of denoting
#'    \code{function(x) x ^ 2}.
#' @param M An \eqn{n\times n} annihilator matrix. If \code{NULL}
#'    (the default), this will be calculated from the \code{mainlm} object
#' @param cluster A logical; should the design matrix X be replaced with an
#'    \eqn{n\times n_c} matrix of ones and zeroes, with a single one in each
#'    row, indicating assignments of the \eqn{n} observations to \eqn{n_c}
#'    clusters using an agglomerative hierarchical clustering algorithm. In
#'    this case, the dimensionality of \eqn{\gamma} is \eqn{n_c} and not
#'    \eqn{p}. Defaults to \code{FALSE}
#' @param varselect Either a character indicating how variable selection should
#'    be conducted, or an integer vector giving indices of columns of the
#'    predictor matrix (\code{\link[stats]{model.matrix}} of \code{mainlm})
#'    to select. The vector must include \code{1L} for the intercept to be
#'    selected. If a character, it must be one of the following:
#' \itemize{
#'  \item \code{"none"}: No variable selection is conducted;
#'  \item \code{"hettest"}: Variable selection is conducted by applying a
#'    heteroskedasticity test with each feature in turn serving as the
#'    `deflator' variable
#'  \item \code{"cv.linear"}: Variable selection is conducted by best subset
#'    selection on the auxiliary linear variance model (linear specification),
#'    using squared-error loss computed under \eqn{K}-fold cross-validation
#'  \item \code{"cv.cluster"}: Variable selection is conducted by best subset
#'    selection on the auxiliary linear variance model (clustering
#'    specification), using squared-error loss computed under \eqn{K}-fold
#'    cross-validation
#'  \item \code{"qgcv.linear"}: Variable selection is conducted by best subset
#'    selection on the auxiliary linear variance model (linear specification),
#'    using squared-error loss computed under quasi-generalised
#'    cross-validation
#'  \item \code{"qgcv.cluster"}: Variable selection is conducted by best subset
#'    selection on the auxiliary linear variance model (clustering
#'    specification), using squared-error loss computed under
#'    quasi-generalised cross-validation
#' }
#' @param nclust A character indicating which elbow method to use to select
#'    the number of clusters (ignored if \code{cluster} is \code{FALSE}).
#'    Alternatively, an integer specifying the number of clusters
#' @param clustering A list object of class \code{"doclust"}. If set to
#'    \code{NULL} (the default), such an object is generated (ignored if
#'    \code{cluster} is \code{FALSE})
#' @param param.init Specifies the initial values of the parameter vector to
#'    use in the Gauss-Newton fitting algorithm. This can either be a function
#'    for generating the initial values from a probability distribution, a
#'    list containing named objects corresponding to the arguments of
#'    \code{\link[base]{seq}} (specifying a sequence of scalar values that
#'    will be passed to \code{\link[base]{expand.grid}}), or a numeric vector
#'    specifying a single initial parameter vector
#' @param maxgridrows An integer indicating the maximum number of initial
#'    values of the parameter vector to try, in case of \code{param.init}
#'    being a function or a list used to generate a grid. Defaults to
#'    \code{20L}.
#' @param nconvstop An integer indicating how many times the quasi-likelihood
#'    estimation algorithm should converge before the grid search across
#'    different initial parameter values is truncated. Defaults to \code{3L}.
#'    If \code{nconvstop >= maxgridrows}, no early stopping rule will be used.
#' @param zerosallowed A logical indicating whether 0 values are acceptable
#'    in the initial values of the parameter vector. Defaults to \code{FALSE}.
#' @param maxitql An integer specifying the maximum number of iterations to
#'    run in the Gauss-Newton algorithm for quasi-likelihood estimation.
#'    Defaults to \code{100L}.
#' @param tolql A double specifying the convergence criterion for the
#'    Gauss-Newton algorithm; defaults to \code{1e-8}. The criterion is applied
#'    to the \code{L_2} norm of the difference between parameter vectors in
#'    successive iterations.
#' @param nestedql A logical indicating whether to use the nested updating step
#'    suggested in \insertCite{Seber03;textual}{skedastic}. Defaults to
#'    \code{FALSE} due to the large computation time required.
#' @param reduce2homosked A logical indicating whether the homoskedastic
#'    error variance estimator \eqn{e'e/(n-p)} should be used if the
#'    variable selection procedure does not select any variables. Defaults to
#'    \code{TRUE}.
#' @param ... Other arguments that can be passed to (non-exported) helper
#'    functions, namely:
#' \itemize{
#'  \item \code{greedy}, a logical passed to the functions implementing best subset
#'    selection, indicating whether or not to use a greedy search rather than
#'    exhaustive search for the best subset. Defaults to \code{FALSE}, but
#'    coerced to \code{TRUE} unconditionally if \eqn{p>9}.
#'  \item \code{distmetric}, a character specifying the distance metric to use in
#'    computing distance for the clustering algorithm. Corresponds to the
#'    \code{method} argument of \code{\link[stats]{dist}} and defaults to
#'    \code{"euclidean"}
#'  \item \code{linkage}, a character specifying the linkage rule to use in
#'    agglomerative hierarchical clustering. Corresponds to the \code{method}
#'    argument of \code{\link[stats]{hclust}} and defaults to
#'    \code{"complete"}
#'  \item \code{alpha}, a double specifying the significance level threshold to use
#'    when applying heteroskedasticity test for the purpose of feature
#'    selection in an ALVM; defaults to \code{0.1}
#'  \item \code{testname}, a character corresponding to the name of a function that
#'    performs a heteroskedasticity test. The function must either be one that
#'    takes a \code{deflator} argument or \code{\link{breusch_pagan}}. Defaults
#'    to \code{evans_king}
#' }
#'
#' @inheritParams breusch_pagan
#' @inheritParams alvm.fit
#'
#' @return An object of class \code{"anlvm.fit"}, containing the following:
#' \itemize{
#'  \item \code{coef.est}, a vector of parameter estimates, \eqn{\hat{\gamma}}
#'  \item \code{var.est}, a vector of estimates \eqn{\hat{\omega}} of the error
#'    variances for all observations
#'  \item \code{method}, either \code{"cluster"} or \code{"functionalform"},
#'    depending on whether \code{cluster} was set to \code{TRUE}
#'  \item \code{ols}, the \code{lm} object corresponding to the original linear
#'    regression model
#'  \item \code{fitinfo}, a list containing three named objects, \code{g} (the
#'    heteroskedastic function), \code{Msq} (the elementwise-square of the
#'    annihilator matrix \eqn{M}), \code{Z} (the design matrix used in the
#'    ANLVM, after feature selection if applicable), and \code{clustering}
#'    (a list object with results of the clustering procedure, if applicable).
#'  \item \code{selectinfo}, a list containing two named objects,
#'    \code{varselect} (the value of the eponymous argument), and
#'    \code{selectedcols} (a numeric vector with column indices of \eqn{X}
#'    that were selected, with \code{1} denoting the intercept column)
#' \item \code{qlinfo}, a list containing nine named objects: \code{converged}
#'    (a logical, indicating whether the Gauss-Newton algorithm converged
#'    for at least one initial value of the parameter vector),
#'    \code{iterations} (the number of Gauss-Newton iterations used to
#'    obtain the parameter estimates returned), \code{Smin} (the minimum
#'    achieved value of the objective function used in the Gauss-Newton
#'    routine), and six arguments passed to the function (\code{nested},
#'    \code{param.init}, \code{maxgridrows}, \code{nconvstop},
#'    \code{maxitql}, and \code{tolql})
#' }
#'
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link{alvm.fit}}, \code{\link{avm.ci}}
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' myanlvm <- anlvm.fit(mtcars_lm, g = function(x) x ^ 2,
#'  varselect = "qgcv.linear")
#'

anlvm.fit <- function(mainlm, g, M = NULL, cluster = FALSE,
        varselect = c("none", "hettest", "cv.linear", "cv.cluster",
            "qgcv.linear", "qgcv.cluster"),
        nclust = c("elbow.swd", "elbow.mwd", "elbow.both"), clustering = NULL,
        param.init = function(q) stats::runif(n = q, min = -5, max = 5),
        maxgridrows = 20L, nconvstop = 3L, zerosallowed = FALSE,
        maxitql = 100L, tolql = 1e-8, nestedql = FALSE,
        reduce2homosked = TRUE, cvoption = c("testsetols", "partitionres"),
        nfolds = 5L, ...) {

  if (!any(is.na(nclust)) && !any(is.numeric(nclust))) {
    nclust <- match.arg(nclust,
                c("elbow.swd", "elbow.mwd", "elbow.both"))
  }
  cvoption <- match.arg(cvoption, c("testsetols", "partitionres"))

  if (is.character(g)) {
    if (g == "sq") {
      g <- function(x) x ^ 2
    } else {
      g <- get(g)
    }
  }

  if (is.character(varselect))
    varselect <- match.arg(varselect,
        c("none", "hettest", "cv.linear", "cv.cluster",
        "qgcv.linear", "qgcv.cluster"))
  if (inherits(mainlm, "lm")) {
    X <- stats::model.matrix(mainlm)
    y <- stats::model.response(stats::model.frame(mainlm))
    e <- stats::resid(mainlm)
    p <- ncol(X)
  } else { # other type of list
    processmainlm(m = mainlm)
  }
  colnames(X) <- NULL
  hasinterceptX <- columnof1s(X)
  esq <- e ^ 2
  n <- length(esq)
  if (is.null(M)) {
    M <- diag(n) - X %*% Rfast::spdinv(crossprod(X)) %*% t(X)
  }
  Msq <- M ^ 2
  novarselected <- FALSE
  dots <- list(...)

  if (is.numeric(varselect)) {
    if (!(1 %in% varselect))
      warning("You have not selected the first column of X.\nThis will result in the omission of an intercept, if present.")
    selectedcols <- varselect
    varselect = "manual"
  } else if (varselect == "none") {
    selectedcols <- 1:p # should always include 1 even after selection
  } else if (varselect == "hettest") {
    selectedcols <- unname(hetvarsel(mainlm, hasintercept = hasinterceptX,
                                     ...)$selectedcols)
    selectedcols <- c(1, selectedcols)
  } else if (varselect == "cv.linear") {
    selectedcols <- unname(varsel.cv.linear(fulldat = list("X" = X, "y" = y,
            "esq" = esq, "Msq" = Msq), nfolds = nfolds, cvopt = cvoption,
            ...)$selectedcols)
  } else if (varselect == "cv.cluster") {
    runsel <- varsel.cv.cluster(fulldat = list("X" = X, "y" = y,
            "esq" = esq, "Msq" = Msq), nclust = nclust,
            nfolds = nfolds, cvopt = cvoption, ...)
    selectedcols <- unname(runsel$selectedcols)
    if (nclust == "foldcv") nclust <- runsel$bestnclust
  } else if (varselect == "qgcv.linear") {
    selectedcols <- unname(varsel.qgcv.linear(mainlm = mainlm,
                                              ...)$selectedcols)
  } else if (varselect == "qgcv.cluster") {
    runsel <- varsel.qgcv.cluster(mainlm = mainlm, nclust = nclust, ...)
    selectedcols <- unname(runsel$selectedcols)
    if (nclust == "foldcv") nclust <- runsel$bestnclust
  }

  if (length(selectedcols) == 1L && selectedcols == 1) {
    novarselected <- TRUE
    # warning("No predictor variables selected for model. Homoskedastic approach used.")
  }

  if (novarselected && reduce2homosked) {
    value <- list("coef.est" = NA_real_,
                  "var.est" = rep(sum(esq) / (n - p), n),
                  "method" = "homoskedastic", "ols" = stats::lm(y ~ 0 + X),
                  "fitinfo" = list("g" = identity,
                                   "Msq" = Msq,
                                   "Z" = matrix(data = 1, nrow = n, ncol = 1),
                                   "clustering" = NULL),
                  "trainerr" = sum((esq - rep(sum(esq) / (n - p), n)) ^ 2) / n,
                  "hyperpar" = list("nclust" = NA_integer_,
                                    "nclustmethod" = NA_character_),
                  "selectinfo" = list("varselect" = varselect,
                                      "selectedcols" = selectedcols),
                  "qlinfo" = list("converged" = TRUE,
                                  "iterations" = NA_integer_,
                                  "Smin" = NA_real_,
                                  "nested" = NA,
                                  "param.init" = NA,
                                  "maxgridrows" = NA_integer_,
                                  "nconvstop" = NA_integer_,
                                  "maxitql" = NA_integer_,
                                  "tolql" = NA_real_))
    class(value) <- "anlvm.fit"
    return(value)
  }

  Xreduced <- X[, selectedcols, drop = FALSE]
  colnames(Xreduced) <- NULL
  hasintercept <- columnof1s(Xreduced)
  if (hasintercept[[1]] & cluster)  {
    # drop intercept so it isn't part of X
    # (clustering shouldn't include intercept)
    Xreduced <- Xreduced[, -hasintercept[[2]], drop = FALSE]
  }

  dots <- list(...)

  if (cluster) {

    if (is.null(clustering)) {
      if (is.numeric(nclust)) {
        thenclust <- nclust
        if (ncol(Xreduced) == 0L) {
          theclustering <- randclust(n = n, nclust = thenclust)
        } else {
          theclustering <- doclust(Xreduced, nclust = thenclust, ...)
        }
      } else {
        if (regexpr("elbow", nclust) != -1) {
          theclustering <- doclust(Xreduced, nclust = nclust, ...)
          thenclust <- length(theclustering$s)
        }
      }
    } else if (!is.null(clustering)) {

      if (inherits(clustering, "doclust")) {
        theclustering <- clustering
        thenclust <- length(theclustering$s)
      } else if (inherits(clustering, "doclustperm")) {
        # theclustering is a list of doclust objects
        # for different variable selection situations
        ww <- which(names(clustering) == paste(selectedcols[-1],
                                               collapse = ""))
        theclustering <- clustering[[ww]]
        thenclust <- length(theclustering$s)
      }
      nclust <- "passed"
    }
    Z <- matrix(data = 0, nrow = n, ncol = thenclust)
    for (j in 1:thenclust) {
      Z[1:n %in% theclustering$s[[j]], j] <- 1
    }
  } else {
    Z <- Xreduced
  }
  q <- ncol(Z)

  # Mean function of ANLVM
  f <- function(gam) {
    gvec <- vapply(1:n, function(k) {
      g(sum(Z[k, ] * gam))
    }, 0)
    drop(Msq %*% gvec)
  }

  # Gradient matrix of mean function of ANLVM
  F. <- function(gam) {
    pracma::jacobian(f = f, x0 = gam)
  }

  # Covariance matrix function of errors of ANLVM
  V <- function(gam) {
    Omega <- diag(vapply(1:n, function(k) {
      g(sum(Z[k, ] * gam))
    }, 0))
    MOmegaM <- M %*% Omega %*% M
    2 * MOmegaM * MOmegaM
  }

  # Standardised sum of squared error function to be minimised
  S <- function(gam) {
    if (is.error(Vinv <- Rfast::spdinv(V(gam)))) {
      NaN
    } else {
      as.numeric(t(esq - f(gam)) %*% Vinv %*% (esq - f(gam)))
    }
  }

  # Updates gam parameter vector once
  updater <- function(gam.a1, gam.ab = NULL) {

    if (is.null(gam.ab)) { # non-nested scheme has only gam.a
      gam.ab <- gam.a1
    }

    if (is.error(V.eval.inv <- Rfast::spdinv(V(gam.a1)))) {
      return(rep(NaN, length(gam.a1)))
    }
    F..eval <- F.(gam.ab)
    f.eval <- f(gam.ab)

    if (is.error(quadform.eval <-
                 Rfast::spdinv(t(F..eval) %*% V.eval.inv %*% F..eval))) {
      return(rep(NaN, length(gam.a1)))
    }

    delta <- quadform.eval %*%
      t(F..eval) %*% V.eval.inv %*% (esq - f.eval)

    gam.ab + delta

  }

  # Performs Gauss-Newton implementation of quasi-likelihood estimation
  quasiopt <- function(param.init, maxiter, convcrit,
                       nested, ...) {

    a <- 0
    convval <- Inf
    convvalprev <- Inf
    gam.anew <- param.init
    gam.aold <- Inf
    while (a < maxiter & (convval > convcrit | convvalprev > convcrit)) {
      a <- a + 1
      convvalprev <- convval
      if (!nested) {
        gam.anew <- updater(gam.anew)
      } else {
        b <- 0
        convval.nested <- Inf
        convval.nestedprev <- Inf
        gam.bold <- gam.anew
        gam.bprev <- Inf

        while (b < maxiter & (convval.nested > convcrit |
                              convval.nestedprev > convcrit)) {
          b <- b + 1
          convval.nestedprev <- convval.nested

          gam.bnew <- updater(gam.anew, gam.bold)

          convval.nested <- sqrt(sum((gam.bnew - gam.bold) ^ 2))

          if (is.na(convval.nested)) {
            convval.nested <- Inf
            break
          }
          gam.bold <- gam.bnew
        }

        if (b >= maxiter | is.infinite(convval.nested))
          warning("At least one instance of nested procedure did not converge")

        gam.anew <- updater(gam.bnew)
      }

      convval <- sqrt(sum((gam.anew - gam.aold) ^ 2))
      if (is.na(convval) || is.nan(convval)) { # results from singular V matrix
        convval <- Inf
        break # abort algorithm and return non-convergence
      }

      gam.aold <- gam.anew
    }

    if (is.na(convval) || is.nan(convval)) {
      didconverge <- FALSE
    } else {
      didconverge <- (convval <= convcrit & convvalprev <= convcrit)
    }
    list("SSE.min" = S(gam.anew), "gam.est" = gam.anew,
         "lastnormdiff" = convval, "numiter.a" = a,
         "converged" = didconverge)
  }

  if (is.numeric(param.init)) { # just use given initial value vector
    if (g(0) == 0 & all(param.init == 0)) {
      stop("Initial values cannot all be 0 with this choice of g function")
    }
    output <- quasiopt(param.init = param.init, maxiter = maxitql,
                       convcrit = tolql, nested = nestedql, ...)
    achievedconvergence <- output$converged
    numiter <- output$numiter.a
    gam.est <- output$gam.est
    Smin <- output$SSE.min
  } else if (is.list(param.init) || is.function(param.init)) {
    # grid search over choice of initial values
    if (is.list(param.init)) {
      if ("by" %in% names(param.init)) {
        allseq <- lapply(1:q, function(j) {
          thisseq <- seq(from = param.init$from, to = param.init$to,
                         by = param.init$by)
          if (!zerosallowed) thisseq[thisseq != 0]
        })
      } else {
        allseq <- lapply(1:q, function(j) {
          thisseq <- seq(from = param.init$from, to = param.init$to,
                         length.out = param.init$length.out)
          if (!zerosallowed) thisseq[thisseq != 0]
        })
      }
      mygrid <- as.matrix(expand.grid(allseq))
      if (g(0) == 0) {
        allzeros <- which(vapply(1:nrow(mygrid), function(i) {
          all(mygrid[i, ] == 0)
        }, TRUE))
        if (length(allzeros) != 0) mygrid <- mygrid[-allzeros, ]
      }
      if (nrow(mygrid) > maxgridrows) {
        mygrid <- mygrid[sample(nrow(mygrid), maxgridrows), ]
      } else { # just randomise order
        mygrid <- mygrid[sample(nrow(mygrid)), ]
      }
    } else if (is.function(param.init)) {
      mygrid <- t(sapply(1:maxgridrows, function(m) {
        param.init(q = q)
      }, simplify = "matrix"))
    }

    alloutput <- vector("list", nrow(mygrid))
    nconv <- 0
    for (m in 1:nrow(mygrid)) {
      alloutput[[m]] <- quasiopt(param.init = mygrid[m, ], maxiter = maxitql,
                            convcrit = tolql, nested = nestedql, ...)
      if (alloutput[[m]]$converged) nconv <- nconv + 1
      if (nconv == nconvstop) break
    }

    allSSEmin <- vapply(1:m, function(mm) {
      alloutput[[mm]]$SSE.min
    }, 0)
    bestsetting <- which.min(allSSEmin)
    if (length(bestsetting) == 0L) {
      numiter <- alloutput[[1]]$numiter.a
      achievedconvergence <- FALSE
      gam.est <- alloutput[[1]]$gam.est
      Smin <- NA_real_
    } else {
      numiter <- alloutput[[bestsetting[1]]]$numiter.a
      achievedconvergence <- alloutput[[bestsetting[1]]]$converged
      gam.est <- alloutput[[bestsetting[1]]]$gam.est
      Smin <- min(allSSEmin)
    }
  }

  value <- list("coef.est" = drop(gam.est),
                "var.est" = vapply(1:n, function(k) {
                  g(sum(Z[k, ] * gam.est))
                }, 0))

  if (any(value$var.est <= 0)) {
    stop("Zero or -ve variance estimates obtained")
  }

  if (!cluster) theclustering <- NULL
  value$method <- ifelse(cluster, "cluster", "functionalform")
  value$ols <- stats::lm(y ~ 0 + X)
  value$fitinfo <- list("g" = g, "Msq" = Msq, "Z" = Z,
                        "clustering" = theclustering)
  value$trainerr <- sum((esq - drop(Msq %*% value$var.est)) ^ 2) / n
  value$hyperpar <- list("nclust" = ifelse(cluster, thenclust, NA_integer_),
                         "nclustmethod" = ifelse(cluster, nclust,
                                                 NA_character_))
  value$selectinfo <- list("varselect" = varselect,
                           "selectedcols" = selectedcols)
  value$qlinfo <- list("converged" = achievedconvergence,
                       "iterations" = numiter,
                       "Smin" = Smin,
                       "nested" = nestedql,
                       "param.init" = param.init,
                       "maxgridrows" = maxgridrows,
                       "nconvstop" = nconvstop,
                       "maxitql" = maxitql,
                       "tolql" = tolql)
  class(value) <- "anlvm.fit"
  return(value)
}
