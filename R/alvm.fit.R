#' Auxiliary Linear Variance Model
#'
#' Fits an Auxiliary Linear Variance Model (ALVM) to estimate the error
#'    variances of a heteroskedastic linear regression model.
#'
#' @details The ALVM model equation is
#'    \deqn{e\circ e = (M \circ M)L \gamma + u},
#'    where \eqn{e} is the Ordinary Least Squares residual vector, \eqn{M} is
#'    the annihilator matrix \eqn{M=I-X(X'X)^{-1}X'}, \eqn{L} is a linear
#'    predictor matrix, \eqn{u} is a random error vector, \eqn{\gamma} is a
#'    \eqn{p}-vector of unknown parameters, and \eqn{\circ} denotes the
#'    Hadamard (elementwise) product. The construction of \eqn{L} depends on
#'    the method used to model or estimate the assumed heteroskedastic
#'    function \eqn{g(\cdot)}, a continuous, differentiable function that is
#'    linear in \eqn{\gamma} and by which the error variances \eqn{\omega_i}
#'    of the main linear model are related to the predictors \eqn{X_{i\cdot}}.
#'    This method has been developed as part of the author's doctoral research
#'    project.
#'
#' Depending on the model used, the estimation method could be
#'    Inequality-Constrained Least Squares or Inequality-Constrained Ridge
#'    Regression. However, these are both special cases of Quadratic
#'    Programming. Therefore, all of the models are fitted using Quadratic
#'    Programming.
#'
#' Several techniques are available for feature selection within the model.
#'    The LASSO-type model handles feature selection via a shrinkage penalty.
#'    For this reason, if the user calls the polynomial model with
#'    \eqn{L_1}-norm penalty, it is not necessary to specify a variable
#'    selection method, since this is handled automatically. Another feature
#'    selection technique is to use a heteroskedasticity test that tests for
#'    heteroskedasticity linked to a particular predictor variable (the
#'    `deflator'). This test can be conducted with each features in turn
#'    serving as the deflator. Those features for which the null hypothesis of
#'    homoskedasticity is rejected at a specified significance level
#'    \code{alpha} are selected. A third feature selection technique is best
#'    subset selection, where the model is fitted with all possible subsets of
#'    features. The models are scored in terms of some metric, and the
#'    best-performing subset of features is selected. The metric could be
#'    squared-error loss computed under \eqn{K}-fold cross-validation or using
#'    quasi-generalised cross-validation. (The \emph{quasi-} prefix refers to
#'    the fact that generalised cross-validation is, properly speaking, only
#'    applicable to a linear fitting method, as defined by
#'    \insertCite{Hastie09;textual}{skedastic}. ALVMs are not linear fitting
#'    methods due to the inequality constraint). Since best subset selection
#'    requires fitting \eqn{2^{p-1}} models (where \eqn{p-1} is the number of
#'    candidate features), it is infeasible for large \eqn{p}. A greedy search
#'    technique can therefore be used as an alternative, where one begins with
#'    a null model and adds the feature that leads to the best improvement in
#'    the metric, stopping when no new feature leads to an improvement.
#'
#' The polynomial and thin-plate spline ALVMs have a penalty hyperparameter
#'    \eqn{\lambda} that must either be specified or tuned. \eqn{K}-fold
#'    cross-validation or quasi-generalised cross-validation can be used for
#'    tuning. The clustering ALVM has a hyperparameter \eqn{n_c}, the number of
#'    clusters into which to group the observations (where error variances
#'    are assumed to be equal within each cluster). \eqn{n_c} can be specified
#'    or tuned. The available tuning methods are an elbow method (using a
#'    sum of within-cluster distances criterion, a maximum
#'    within-cluster distance criterion, or a combination of the two) and
#'    \eqn{K}-fold cross-validation.
#'
#' @param model A character corresponding to the type of ALVM to be fitted:
#'    \code{"cluster"} for the clustering ALVM, \code{"spline"} for the
#'    thin-plate spline ALVM, \code{"linear"} for the linear ALVM,
#'    \code{"polynomial"} for the penalised polynomial ALVM, \code{"basic"} for
#'    the basic or na{\" i}ve ALVM, and \code{"homoskedastic"} for the
#'    homoskedastic error variance estimator, \eqn{e'e/(n-p)}.
#' @param lambda Either a double of length 1 indicating the value of the
#'    penalty hyperparameter \eqn{\lambda}, or a character specifying the
#'    tuning method for choosing \eqn{\lambda}: \code{"foldcv"} for
#'    \eqn{K}-fold cross-validation (the default) or \code{"qgcv"} for
#'    quasi-generalised cross-validation. This argument is ignored if
#'    \code{model} is neither \code{"polynomial"} nor \code{"spline"}.
#' @param nclust Either an integer of length 1 indicating the value of the
#'    number of clusters \eqn{n_c}, or a character specifying the tuning method
#'    for choosing \eqn{n_c}: \code{"elbow.swd"} for the elbow method using a
#'    sum of within-cluster distances criterion,  \code{"elbow.mwd"} for the
#'    elbow method using a maximum within-cluster distances criterion,
#'    \code{"elbow.both"} for rounded average of the results of the previous
#'    two, and \code{"foldcv"} for \eqn{K}-fold cross-validation. This argument
#'    is ignored if \code{model} is not \code{"cluster"}.
#' @param polypen A character, either \code{"L2"} or \code{"L1"}, indicating
#'    whether an \eqn{L_2} norm penalty (ridge regression) or
#'    \eqn{L_1} norm penalty (LASSO) should be used with the polynomial model.
#'    This argument is ignored if \code{model} is not \code{"polynomial"}.
#' @param solver A character, indicating which Quadratic Programming solver
#'    function to use to estimate \eqn{\gamma}. The options are
#'    \code{"quadprog"}, corresponding to
#'    \code{\link[quadprog]{solve.QP.compact}} from the \pkg{quadprog} package;
#'    package; \code{"quadprogXT"}, corresponding to
#'    \code{\link[quadprogXT]{buildQP}} from the \pkg{quadprogXT} package;
#'    \code{"roi"}, corresponding to the \code{qpoases} solver implemented in
#'    \code{\link[ROI]{ROI_solve}} from the \pkg{ROI} package with
#'    \pkg{ROI.plugin.qpoases} add-on; and \code{"osqp"}, corresponding to
#'    \code{\link[osqp]{solve_osqp}} from the \pkg{osqp} package.
#'    Alternatively, the user can specify \code{"auto"} (the default), in which
#'    case the function will select the solver that seems to work best for the
#'    chosen model.
#' @param constol A double corresponding to the boundary value for the
#'    constraint on error variances. Of course, the error variances must be
#'    non-negative, but setting the constraint boundary to 0 can result in
#'    zero estimates that then result in infinite weights for Feasible
#'    Weighted Least Squares. The boundary value should thus be positive, but
#'    small enough not to bias estimation of very small variances. Defaults to
#'    \code{1e-10}.
#' @param tsk An integer corresponding to the basis dimension \code{k} to be
#'    passed to the \code{[mgcv]{s}} function for fitting of a thin-plate
#'    spline ALVM; see \code{[mgcv]{choose.k}} for more details about
#'    choosing the parameter and defaults. Ignored if \code{model} is not
#'    \code{"spline"}.
#' @param tsm An integer corresponding to the order \code{m} of the penalty
#'    to be passed to the \code{[mgcv]{s}} function for fitting of a thin-plate
#'    spline ALVM. If left as the default (\code{NULL}), it will be set to
#'    2, corresponding to 2nd derivative penalties for a cubic spline.
#' @param cvoption A character, either \code{"testsetols"} or
#'    \code{"partitionres"}, indicating how to obtain the observed response
#'    values for each test fold when performing \eqn{K}-fold cross-validation
#'    on an ALVM. The default technique, \code{"testsetols"}, entails fitting
#'    a linear regression model to the test fold of observations from the
#'    original response vector \eqn{y} and predictor matrix \eqn{X}. The
#'    squared residuals from this regression are the observed
#'    responses that are predicted from the trained model to compute the
#'    cross-validated squared error loss function. Under the other technique,
#'    \code{"partitionres"}, the squared residuals from the full
#'    linear regression model are partitioned into training and test folds and
#'    the squared residuals in the test fold are the observed responses that
#'    are predicted for computation of the cross-validated loss.
#' @param nfolds An integer specifying the number of folds \eqn{K} to use for
#'    cross-validation, if the \eqn{\lambda} and/or \eqn{n_c} hyperparameters
#'    are to be tuned using cross-validation. Defaults to \code{5L}. One must
#'    ensure that each test fold contains at least \eqn{p+1} observations if
#'    the \code{"testsetols"} technique is used with cross-validation, so that
#'    there are enough degrees of freedom to fit a linear model to the test
#'    fold.
#' @param d An integer specifying the degree of polynomial to use in the
#'    penalised polynomial ALVM; defaults to \code{2L}. Ignored if
#'    \code{model} is other than \code{"polynomial"}. Setting \code{d} to
#'    \code{1L} is not identical to setting \code{model} to \code{"linear"},
#'    because the linear ALVM does not have a penalty term in the objective
#'    function.
#' @param ... Other arguments that can be passed to (non-exported) helper
#'    functions, namely:
#' \itemize{
#'  \item \code{greedy}, a logical passed to the functions implementing best subset
#'    selection, indicating whether or not to use a greedy search rather than
#'    exhaustive search for the best subset. Defaults to \code{FALSE}, but
#'    coerced to \code{TRUE} unconditionally if \eqn{p>9}.
#' \item \code{distmetric}, a character specifying the distance metric to use in
#'    computing distance for the clustering algorithm. Corresponds to the
#'    \code{method} argument of \code{\link[stats]{dist}} and defaults to
#'    \code{"euclidean"}
#' \item \code{linkage}, a character specifying the linkage rule to use in
#'    agglomerative hierarchical clustering. Corresponds to the \code{method}
#'    argument of \code{\link[stats]{hclust}} and defaults to
#'    \code{"complete"}
#' \item \code{nclust2search}, an integer vector specifying the values of
#'    \eqn{n_c} to try when tuning \eqn{n_c} by cross-validation. Defaults to
#'    \code{1L:20L}
#' \item \code{alpha}, a double specifying the significance level threshold to
#'    use when applying heteroskedasticity test for the purpose of feature
#'    selection in an ALVM; defaults to \code{0.1}
#' \item \code{testname}, a character corresponding to the name of a function
#'    that performs a heteroskedasticity test. The function must either be one
#'    that takes a \code{deflator} argument or \code{\link{breusch_pagan}}.
#'    Defaults to \code{evans_king}
#' }
#'
#' @inheritParams breusch_pagan
#' @inheritParams anlvm.fit
#'
#' @importFrom mgcv s
#' @import ROI.plugin.qpoases
#'
#' @return An object of class \code{"alvm.fit"}, containing the following:
#' \itemize{
#'  \item \code{coef.est}, a vector of parameter estimates, \eqn{\hat{\gamma}}
#'  \item \code{var.est}, a vector of estimates \eqn{\hat{\omega}} of the error
#'    variances for all observations
#'  \item \code{method}, a character corresponding to the \code{model} argument
#'  \item \code{ols}, the \code{lm} object corresponding to the original linear
#'    regression model
#'  \item \code{fitinfo}, a list containing four named objects: \code{Msq} (the
#'    elementwise-square of the annihilator matrix \eqn{M}), \code{L} (the
#'    linear predictor matrix \eqn{L}), \code{clustering} (a list object
#'    with results of the clustering procedure), and \code{gam.object}, an
#'    object of class \code{"gam"} (see \code{\link[mgcv]{gamObject}}). The
#'    last two are set to \code{NA} unless the clustering ALVM or thin-plate
#'    spline ALVM is used, respectively
#'  \item \code{hyperpar}, a named list of hyperparameter values,
#'    \code{lambda}, \code{nclust}, \code{tsk}, and \code{d}, and tuning
#'    methods, \code{lambdamethod} and \code{nclustmethod}. Values
#'    corresponding to unused hyperparameters are set to \code{NA}.
#'  \item \code{selectinfo}, a list containing two named objects,
#'    \code{varselect} (the value of the eponymous argument), and
#'    \code{selectedcols} (a numeric vector with column indices of \eqn{X}
#'    that were selected, with \code{1} denoting the intercept column)
#'  \item \code{pentype}, a character corresponding to the \code{polypen}
#'    argument
#'  \item \code{solver}, a character corresponding to the \code{solver}
#'    argument (or specifying the QP solver actually used, if \code{solver}
#'    was set to \code{"auto"})
#'  \item \code{constol}, a double corresponding to the \code{constol} argument
#' }
#'
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link{alvm.fit}}, \code{\link{avm.ci}}
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' myalvm <- alvm.fit(mtcars_lm, model = "polynomial", polypen = "L1")
#'

alvm.fit <- function(mainlm, M = NULL,
      model = c("cluster", "spline", "linear",
        "polynomial", "basic", "homoskedastic"),
    varselect = c("none", "hettest", "cv.linear", "cv.cluster",
        "qgcv.linear", "qgcv.cluster"),
      lambda = c("foldcv", "qgcv"),
        nclust = c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"),
        clustering = NULL, polypen = c("L2", "L1"), d = 2L,
        solver = c("auto", "quadprog",
            "quadprogXT", "roi", "osqp"),
        tsk = NULL, tsm = NULL, constol = 1e-10,
        cvoption = c("testsetols", "partitionres"),
        nfolds = 5L, reduce2homosked = TRUE, ...) {

  model <- match.arg(model, c("cluster", "spline", "linear",
                              "polynomial", "basic", "homoskedastic"))
  solver <- match.arg(solver, c("auto", "quadprog",
                                "quadprogXT", "roi", "osqp"))
  polypen <- match.arg(polypen, c("L2", "L1"))
  cvoption <- match.arg(cvoption, c("testsetols", "partitionres"))
  if (!any(is.na(lambda)) && !any(is.numeric(lambda))) {
    lambda <- match.arg(lambda, c("foldcv", "qgcv"))
  }
  if (!any(is.na(nclust)) && !any(is.numeric(nclust))) {
    nclust <- match.arg(nclust,
        c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"))
  }

  if (solver == "auto") {
    if (model %in% c("cluster", "linear", "basic")) {
      thesolver <- "quadprog"
    } else if ((model == "spline") ||
               (model == "polynomial" && polypen == "L2")) {
      thesolver <- "osqp"
    } else if (model == "polynomial" && polypen == "L1") {
      thesolver <- "roi"
    } else if (model == "homoskedastic") {
      thesolver <- "none"
    }
  } else {
    thesolver <- solver
  }
  if (is.numeric(lambda)) {
    thelambda <- lambda
  } else if (is.null(lambda)) {
    thelambda <- 0
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
            "esq" = esq, "Msq" = Msq), nfolds = nfolds,
            cvopt = cvoption, ...)$selectedcols)
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

  if (model == "homoskedastic" || (novarselected && reduce2homosked)) {
    value <- list("coef.est" = NA_real_,
                  "var.est" = rep(sum(esq) / (n - p), n),
                  "method" = "homoskedastic", "ols" = stats::lm(y ~ 0 + X),
                  "fitinfo" = list("Msq" = Msq, "L" = diag(n),
                      "clustering" = NA, "gam.object" = NA),
                  "hyperpar" = list("lambda" = NA_real_,
                                    "nclust" = NA_integer_,
                        "d" = NA_integer_, "tsk" = NA_integer_,
                        "lambdamethod" = lambda, "nclustmethod" = nclust),
                  "selectinfo" = list("varselect" = varselect,
                                      "selectedcols" = selectedcols),
                  "pentype" = polypen,
                  "solver" = thesolver,
                  "constol" = constol)
    class(value) <- "alvm.fit"
    return(value)
  }

  Xreduced <- X[, selectedcols, drop = FALSE]
  colnames(Xreduced) <- NULL
  q <- ncol(Xreduced)
  hasintercept <- columnof1s(Xreduced)
  if (hasintercept[[1]] && (model %in% c("spline", "cluster")))  {
    # drop intercept so it isn't part of X
    # (spline functions introduce intercept separately;
    # clustering shouldn't include intercept)
    Xreduced <- Xreduced[, -hasintercept[[2]], drop = FALSE]
  }

  datvars <- list("X" = X, "Xreduced" = Xreduced, "Msq" = Msq,
                  "y" = y, "esq" = esq)
  dots <- list(...)

  if (model == "basic") {
    L <- A <- diag(n)
    P <- matrix(data = 0, nrow = n, ncol = n)
    thepentype <- "L2"
    thelambda <- 0
  } else if (model == "linear") {
    thed <- ifelse("d" %in% names(dots), dots$d, 1L)
    if (thed == 1L) {
      L <- A <- Xreduced
    } else {
      L <- A <- makepolydesign(Xreduced, thed)
    }
    P <- matrix(data = 0, nrow = ncol(L), ncol = ncol(L))
    thelambda <- 0
    thepentype <- "L2"
  } else if (model == "cluster") {
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
        } else if (nclust == "foldcv") {
          thenclust <- tune.nclust(fulldat = datvars,
            solver = thesolver, cvopt = cvoption, nfolds = nfolds, ...)
          theclustering <- doclust(Xreduced, nclust = thenclust, ...)
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
    thepentype <- "L2"
    thelambda <- 0
    L <- matrix(data = 0, nrow = n, ncol = thenclust)
    for (j in 1:thenclust) {
      L[1:n %in% theclustering$s[[j]], j] <- 1
    }
    A <- diag(thenclust)
    P <- matrix(data = 0, nrow = thenclust, ncol = thenclust)
  } else if (model == "spline") {
    thindat <- data.frame("esq" = esq, "X" = Xreduced)
    if (!is.null(tsk) && !is.null(tsm)) {
      argbits <- paste0(", k = ", tsk, ", m = ", tsm)
    } else if (!is.null(tsk) && is.null(tsm)) {
      argbits <- paste0(", k = ", tsk)
    } else if (is.null(tsk) && !is.null(tsm)) {
      argbits <- paste0(", m = ", tsm)
    } else {
      argbits <- character(0)
    }
    if (ncol(Xreduced) > 1) {
      thinformula <- stats::as.formula(paste0("esq ~ s(",
            paste0("X.", 1:ncol(Xreduced), collapse = ","),
            ", bs = \"tp\"", argbits, ")"))
    } else {
      thinformula <- stats::as.formula(paste0("esq ~ s(X, bs = \"tp\"",
            argbits, ")"))
    }
    gam.thin <- mgcv::gam(thinformula, data = thindat, ...)
    Dmat <- cbind(1,
              mgcv::PredictMat(gam.thin$smooth[[1]], data = thindat))
    L <- A <- Rfast::spdinv(Msq) %*% Dmat
    sobj <- eval(parse(text = as.character(thinformula)[3]))
    P <- mgcv::smooth.construct.tp.smooth.spec(sobj,
                  data = thindat, knots = NULL)$S[[1]]
    thepentype <- "L2"
    if (lambda == "foldcv") {
      thelambda <- tune.lambda.cv(fulldat = datvars, method = "spline",
            solver = thesolver, pentype = "L2", tsk = tsk, tsm = tsm,
            cvopt = cvoption, Dmat = Dmat, nfolds = nfolds, ...)
    } else if (lambda == "qgcv") {
      thelambda <- tune.lambda.qgcv(mats = list("D" = Dmat,
          "Dcross" = crossprod(Dmat), "esq" = esq, "A" = A,
          "P" = P, "bvec" = rep(constol, n), "cvec" = drop(t(Dmat) %*% esq)),
          solver = thesolver, ...)
    }
  } else if (model == "polynomial") {
    thed <- d
    L <- A <- makepolydesign(Xreduced, thed)
    thepentype <- polypen
    Dmat <- Msq %*% L
    if (polypen == "L1") {
      P <- rep(1, 2 * ncol(L))
      if (hasintercept[[1]]) P[c(1, ncol(L) + 1)] <- 0
    } else if (polypen == "L2") {
      P <- diag(ncol(L))
      if (hasintercept[[1]]) P[1, 1] <- 0
    }
    if (lambda == "foldcv") {
      thelambda <- tune.lambda.cv(fulldat = datvars, method = "polynomial",
        solver = thesolver, pentype = polypen, cvopt = cvoption,
        Dmat = Dmat, nfolds = nfolds, d = thed, ...)
    } else if (lambda == "qgcv") {
      if (polypen == "L2") {
        thelambda <- tune.lambda.qgcv(mats = list("D" = Dmat,
          "Dcross" = crossprod(Dmat), "esq" = esq, "A" = A,
          "P" = P, "bvec" = rep(constol, n), "cvec" = drop(t(Dmat) %*% esq)),
          solver = thesolver, ...)
      } else if (polypen == "L1") {
        thelambda <- tune.lambda.qgcv.lasso(mats = list("D" = Dmat,
          "Dcross" = crossprod(Dmat), "esq" = esq, "A" = A,
          "P" = P, "bvec" = rep(constol, n),
          "cvec" = drop(t(Dmat) %*% esq)), ...)
      }
    }
  }

  if (!exists("Dmat", inherits = FALSE)) Dmat <- Msq %*% L
  Dcross <- crossprod(Dmat)
  cvec <- drop(t(Dmat) %*% esq)
  if (model == "cluster") {
    bveclength <- thenclust
  } else {
    bveclength <- n
  }

  qpsol <- qpest(Dcross = Dcross, P = P, A = A,
                 cvec = cvec, bvec = rep(constol, bveclength),
                 lambda = thelambda, solver = thesolver, ...)

  value <- list("coef.est" = qpsol$coef.est,
                "var.est" = drop(L %*% qpsol$coef.est))

  if (any(value$var.est <= 0)) {
    negzeroest <- TRUE
    warning("Zero or -ve variance estimates obtained; trying other solvers...")
    solverpool <- setdiff(c("quadprog", "osqp", "roi"),
                          thesolver)
    for (slvr in solverpool) {
      qpsol <- qpest(Dcross = Dcross, P = P, A = A,
                     cvec = cvec, bvec = rep(constol, bveclength),
                     lambda = thelambda, solver = slvr, ...)
      value <- list("coef.est" = qpsol$coef.est,
                    "var.est" = drop(L %*% qpsol$coef.est))
      if (all(value$var.est > 0)) {
        negzeroest <- FALSE
        thesolver <- slvr
        break
        message(paste0(thesolver, " solver successful."))
        theusecov <- "no"
      }
    }
    if (negzeroest)
      warning("Zero or -ve variance estimates obtained; Trying different solvers did not rectify.")
  }

  if (model != "cluster") theclustering <- NA
  value$method <- model
  value$ols <- stats::lm(y ~ 0 + X)
  value$fitinfo <- list("Msq" = Msq, "L" = L, "clustering" = theclustering)
  value$fitinfo$gam.object <- if (model == "spline") {
    gam.thin
  } else {
    NA
  }
  value$hyperpar$lambda <- if (model %in% c("polynomial", "spline")) {
    thelambda
  } else {
    NA_real_
  }
  value$hyperpar$nclust <- ifelse(model == "cluster", thenclust, NA_integer_)
  value$hyperpar$d <- if (model == "polynomial") {
    thed
  } else {
    NA_integer_
  }
  value$hyperpar$tsk <- if (model == "spline") {
    ifelse(is.null(tsk), length(value$coef.est), tsk)
  } else {
    NA_integer_
  }

  value$hyperpar$lambdamethod <- lambda
  value$hyperpar$nclustmethod <- ifelse(model == "cluster",
                                        nclust, NA_character_)
  value$selectinfo <- list("varselect" = varselect,
                           "selectedcols" = selectedcols)
  value$pentype <- polypen
  value$solver <- thesolver
  value$constol <- constol
  class(value) <- "alvm.fit"
  return(value)
}
