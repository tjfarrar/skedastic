varsel.qgcv.linear <- function(mainlm, greedy = FALSE, ...) {

  processmainlm(m = mainlm)
  esq <- e ^ 2
  n <- length(esq)
  if (!greedy && p > 9) {
    warning("greedy changed to TRUE since 2 ^ (p-1) > 100")
    greedy <- TRUE
  }
  allsubsets <- as.matrix(expand.grid(replicate(p - 1, c(FALSE, TRUE),
      simplify = FALSE)))
  allsubsets <- cbind(TRUE, allsubsets)

  if (greedy) { # greedy search for best subset
    bestObj.smallermodel <- Inf
    bestcovariates <- 1
    for (pprime in 1:p) { # number of vars (incl. intercept) in model
      thissubsets <- allsubsets[rowSums(allsubsets) == pprime, ,
                                drop = FALSE]
      has1inrightplaces <- vapply(1:nrow(thissubsets), function(l) {
        all(thissubsets[l, bestcovariates] == 1)
      }, TRUE)
      thissubsets <- thissubsets[has1inrightplaces, , drop = FALSE]

      qgcvloss <- vapply(1:nrow(thissubsets), function(l) {
        fit <- alvm.fit(mainlm = mainlm, varselect = which(thissubsets[l, ]),
                        model = "linear", reduce2homosked = FALSE, ...)
        D <- fit$fitinfo$Msq %*% fit$fitinfo$L
        HD <- D %*% Rfast::spdinv(crossprod(D)) %*% t(D)
        trHD <- sum(diag(HD))
        esqpred <- drop(fit$fitinfo$Msq %*% fit$var.est)
        sum(((esq - esqpred) / (1 - trHD / n)) ^ 2) / n
      }, 0)

      bestObj <- min(qgcvloss)
      if (bestObj > bestObj.smallermodel) {
        break # best model of this size is worse than best model of smaller size
      } else {
        bestObj.smallermodel <- bestObj
        bestcovariates <- which(thissubsets[which.min(qgcvloss), ])
      }
    }
    result <- list("selectedcols" = bestcovariates)
  } else if (!greedy) { # exhaustive search for best subset

    qgcvloss <- vapply(1:nrow(allsubsets), function(l) {
      fit <- alvm.fit(mainlm = mainlm, varselect = which(allsubsets[l, ]),
                model = "linear", reduce2homosked = FALSE, ...)
      D <- fit$fitinfo$Msq %*% fit$fitinfo$L
      HD <- D %*% Rfast::spdinv(crossprod(D)) %*% t(D)
      trHD <- sum(diag(HD))
      esqpred <- drop(fit$fitinfo$Msq %*% fit$var.est)
      sum(((esq - esqpred) / (1 - trHD / n)) ^ 2) / n
    }, 0)

    bestl <- which.min(qgcvloss)
    result <- list("selectedcols" = which(allsubsets[bestl, ]))
  }
  result
}

varsel.qgcv.cluster <- function(mainlm,
    nclust = c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"),
    greedy = FALSE, ...) {

  if (!any(is.na(nclust)) && !any(is.numeric(nclust))) {
    nclust <- match.arg(nclust,
        c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"))
    if (nclust == "foldcv") nclust <- "elbow.swd"
  }

  processmainlm(m = mainlm)
  esq <- e ^ 2
  n <- length(esq)
  if (!greedy && p > 9) {
    warning("greedy changed to TRUE since 2 ^ (p-1) > 100")
    greedy <- TRUE
  }
  allsubsets <- as.matrix(expand.grid(replicate(p - 1, c(FALSE, TRUE),
        simplify = FALSE)))
  allsubsets <- cbind(TRUE, allsubsets)
  allsubsets <- allsubsets[-1, , drop = FALSE]

  nullmodelval <- vapply(1, function(r) {
    fit <- alvm.fit(mainlm = mainlm, varselect = "none",
                    model = "cluster", nclust = 1L,
                    reduce2homosked = FALSE, ...)
    D <- fit$fitinfo$Msq %*% fit$fitinfo$L
    HD <- D %*% Rfast::spdinv(crossprod(D)) %*% t(D)
    trHD <- sum(diag(HD))
    esqpred <- drop(fit$fitinfo$Msq %*% fit$var.est)
    sum(((esq - esqpred) / (1 - trHD / n)) ^ 2) / n
  }, 0)

  if (greedy) { # greedy search for best subset
    bestObj.smallermodel <- Inf
    bestcovariates <- 1
    bestnclust <- 1L
    for (pprime in 2:p) { # number of vars (incl. intercept) in model
      thissubsets <- allsubsets[rowSums(allsubsets) == pprime, ,
                                drop = FALSE]
      has1inrightplaces <- vapply(1:nrow(thissubsets), function(l) {
        all(thissubsets[l, bestcovariates] == 1)
      }, TRUE)
      thissubsets <- thissubsets[has1inrightplaces, , drop = FALSE]

      qgcvloss <- lapply(1:nrow(thissubsets), function(l) {
        fit <- alvm.fit(mainlm = mainlm, varselect = which(thissubsets[l, ]),
                  model = "cluster", nclust = nclust,
                  reduce2homosked = FALSE, ...)
        D <- fit$fitinfo$Msq %*% fit$fitinfo$L
        HD <- D %*% Rfast::spdinv(crossprod(D)) %*% t(D)
        trHD <- sum(diag(HD))
        esqpred <- drop(fit$fitinfo$Msq %*% fit$var.est)
        list("loss" = sum(((esq - esqpred) / (1 - trHD / n)) ^ 2) / n,
             "bestnclust" = fit$hyperpar$nclust)
      })
      allObj <- vapply(qgcvloss, function(x) x$loss, 0)
      allnclust <- vapply(qgcvloss, function(x) x$bestnclust, 0L)

      minsubs <- which.min(allObj)
      bestObj <- allObj[minsubs]
      if (bestObj > bestObj.smallermodel) {
        break # best model of this size is worse than best model of smaller size
      } else {
        bestObj.smallermodel <- bestObj
        bestcovariates <- which(thissubsets[minsubs, ])
        bestnclust <- allnclust[minsubs]
      }
    }
    result <- list("selectedcols" = bestcovariates,
                   "bestnclust" = bestnclust)
  } else if (!greedy) { # exhaustive search for best subset

    qgcvloss <- lapply(1:nrow(allsubsets), function(l) {
      fit <- alvm.fit(mainlm = mainlm, varselect = which(allsubsets[l, ]),
                      model = "cluster", nclust = nclust,
                      reduce2homosked = FALSE, ...)
      D <- fit$fitinfo$Msq %*% fit$fitinfo$L
      HD <- D %*% Rfast::spdinv(crossprod(D)) %*% t(D)
      trHD <- sum(diag(HD))
      esqpred <- drop(fit$fitinfo$Msq %*% fit$var.est)
      list("loss" = sum(((esq - esqpred) / (1 - trHD / n)) ^ 2) / n,
           "bestnclust" = fit$hyperpar$nclust)
    })
    allObj <- vapply(qgcvloss, function(x) x$loss, 0)
    allnclust <- vapply(qgcvloss, function(x) x$bestnclust, 0L)

    bestl <- which.min(allObj)
    result <- list("selectedcols" = which(allsubsets[bestl, ]),
                   "bestnclust" = allnclust[bestl])

    if (nullmodelval <= allObj[bestl]) { # null model is best
      result <- list("selectedcols" = 1L,
                     "bestnclust" = 1L)
    }
  }
  result
}

# Implemented for linear method only
varsel.cv.linear <- function(fulldat, nfolds = 5L,
              cvopt = "testsetols", greedy = FALSE, ...) {

  p <- ncol(fulldat$X)
  if (!greedy && p > 9) {
    warning("greedy changed to TRUE since 2 ^ (p-1) > 100")
    greedy <- TRUE
  }
  allsubsets <- as.matrix(expand.grid(replicate(p - 1, c(FALSE, TRUE),
                          simplify = FALSE)))
  allsubsets <- cbind(TRUE, allsubsets)
  cvind <- makefolds(fulldat$esq, K = nfolds)

  if (greedy) { # greedy search for best subset (assumes independence of covariates)
    bestObj.smallermodel <- Inf
    bestcovariates <- 1
    for (pprime in 1:p) { # number of vars (incl. intercept) in model
      # Get all candidate variable sets of this size,
      # with already-selected variables strictly included
      thissubsets <- allsubsets[rowSums(allsubsets) == pprime, ,
                                drop = FALSE]
      has1inrightplaces <- vapply(1:nrow(thissubsets), function(l) {
        all(thissubsets[l, bestcovariates] == 1)
      }, TRUE)
      thissubsets <- thissubsets[has1inrightplaces, , drop = FALSE]

      trainlist <- lapply(1:nrow(thissubsets), function(l) {
        fulldat$Xreduced <- fulldat$X[, thissubsets[l, ], drop = FALSE]
        traincalc(fulldat = fulldat, ind = cvind$train, method = "linear",
                  pentype = "L2", ...)
      })
      testlist <- lapply(1:nrow(thissubsets), function(l) {
        fulldat$Xreduced <- fulldat$X[, thissubsets[l, ], drop = FALSE]
        testcalc(fulldat = fulldat, ind = cvind$test,
                 method = "linear", pentype = "L2", cvoption = cvopt, ...)
      })

      myObjFun <- vapply(1:nrow(thissubsets), function(l) {
        Dm <- if (cvopt == "testsetols") {
          NULL
        } else {
          fulldat$Msq %*% fulldat$X[, thissubsets[l, ], drop = FALSE]
        }
        CVObjFun.varsel(train = trainlist[[l]],
            test = testlist[[l]], method = "linear", cvoption = cvopt,
            ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, ...)
      }, 0)
      bestObj <- min(myObjFun)
      if (bestObj > bestObj.smallermodel) {
        break # best model of this size is worse than best model of smaller size
      } else {
        bestObj.smallermodel <- bestObj
        bestcovariates <- which(thissubsets[which.min(myObjFun), ])
      }
    }
    result <- list("selectedcols" = bestcovariates)
  } else if (!greedy) { # exhaustive search for best subset

    trainlist <- lapply(1:nrow(allsubsets), function(l) {
      fulldat$Xreduced <- fulldat$X[, allsubsets[l, ], drop = FALSE]
      traincalc(fulldat = fulldat, ind = cvind$train, method = "linear",
                pentype = "L2", ...)
    })
    testlist <- lapply(1:nrow(allsubsets), function(l) {
      fulldat$Xreduced <- fulldat$X[, allsubsets[l, ], drop = FALSE]
      testcalc(fulldat = fulldat, ind = cvind$test,
               method = "linear", pentype = "L2", cvoption = cvopt, ...)
    })

    myObjFun <- vapply(1:nrow(allsubsets), function(l) {
      Dm <- if (cvopt == "testsetols") {
        NULL
      } else {
        fulldat$Msq %*% fulldat$X[, allsubsets[l, ], drop = FALSE]
      }
      CVObjFun.varsel(train = trainlist[[l]],
                      test = testlist[[l]], method = "linear", cvoption = cvopt,
                      ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, ...)
    }, 0)

    bestl <- which.min(myObjFun)
    result <- list("selectedcols" = which(allsubsets[bestl, ]))
  }
  result
}

varsel.cv.cluster <- function(fulldat,
      nclust = c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"),
      nfolds = 5L, cvopt = "testsetols", greedy = FALSE, ...) {

  if (!any(is.na(nclust)) && !any(is.numeric(nclust))) {
    nclust <- match.arg(nclust,
      c("elbow.swd", "elbow.mwd", "elbow.both", "foldcv"))
    if (nclust == "foldcv") nclust <- "elbow.swd"
  }

  if (nclust == "foldcv") {
    nclust2pass <- NULL
  } else {
    nclust2pass <- nclust
  }
  p <- ncol(fulldat$X)

  if (!greedy && p > 9) {
    warning("greedy changed to TRUE since 2 ^ (p-1) > 100")
    greedy <- TRUE
  }
  allsubsets <- as.matrix(expand.grid(replicate(p - 1, c(FALSE, TRUE),
      simplify = FALSE)))
  allsubsets <- cbind(TRUE, allsubsets)
  allsubsets <- allsubsets[-1, , drop = FALSE]
  cvind <- makefolds(fulldat$esq, K = nfolds)

  # NULL MODEL (only 1 cluster)
  Dm <- if (cvopt == "testsetols") {
    NULL
  } else {
    fulldat$Msq %*% matrix(data = 1, nrow = nrow(fulldat$Msq),
                           ncol = 1)
  }
  train0 <- traincalc0(fulldat = fulldat, ind = cvind$train)
  test0 <- testcalc0(fulldat = fulldat, ind = cvind$test,
                     cvoption = cvopt)
  nullmodelval <- CVObjFun.varsel(train = train0,
      test = test0, method = "cluster", cvoption = cvopt,
      ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, nclust = 1L, ...)

  if (greedy) { # greedy search for best subset (assumes independence of covariates)
    for (pprime in 2:p) { # number of vars in model
      # (THERE IS NO INTERCEPT IN CLUSTER MODEL)
      # Get all candidate variable sets of this size,
      # with already-selected variables strictly included
      bestcovariates <- 1
      bestnclust <- 1L
      bestObj.smallermodel <- nullmodelval
      thissubsets <- allsubsets[rowSums(allsubsets) == pprime, ,
                                drop = FALSE]
      has1inrightplaces <- vapply(1:nrow(thissubsets), function(l) {
        all(thissubsets[l, bestcovariates] == 1)
      }, TRUE)
      thissubsets <- thissubsets[has1inrightplaces, , drop = FALSE]

      cvlist <- lapply(1:nrow(thissubsets), function(l) {
        fulldat$Xreduced <- fulldat$X[, thissubsets[l, ], drop = FALSE]
        train <- traincalc(fulldat = fulldat, ind = cvind$train,
                    method = "cluster", pentype = "L2",
                    nclustmetric = nclust2pass, ...)
        test <- testcalc(fulldat = fulldat, ind = cvind$test,
                    method = "cluster", pentype = "L2", cvoption = cvopt,
                    train = train, ...)
        list("train" = train, "test" = test)
      })

      myObjFun <- lapply(1:nrow(thissubsets), function(l) {
        if (nclust == "foldcv") {
          Dm <- NULL
          maxnclust <- length(cvlist[[l]]$train[[1]]$allclust)
          allCVObjFunvals <- rep(Inf, maxnclust)
          for (r in 2:maxnclust) {
            allCVObjFunvals[r] <- CVObjFun.varsel(train = cvlist[[l]]$train,
             test = cvlist[[l]]$test, method = "cluster", cvoption = cvopt,
             ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, nclust = r, ...)
            if (r >= 4) {
              if (allCVObjFunvals[r] > allCVObjFunvals[r - 1] &&
                  allCVObjFunvals[r] > allCVObjFunvals[r - 2] &&
                  allCVObjFunvals[r] > allCVObjFunvals[r - 3]) break
            }
          }
          list("minfunval" = min(allCVObjFunvals),
               "bestnclust" = which.min(allCVObjFunvals))
        } else {
          Dm <- if (cvopt == "testsetols") {
            NULL
          } else {
            clustering <- doclust(X = fulldat$X[, thissubsets[l, ],
              drop = FALSE], nclust = nclust, ...)
            fulldat$L <- matrix(data = 0, nrow = nrow(fulldat$Msq),
                                ncol = length(clustering$s))
            for (i in 1:nrow(fulldat$Msq)) {
              fulldat$L[i, clustering$C[i]] <- 1
            }
            fulldat$Msq %*% fulldat$L
          }
          funval <- CVObjFun.varsel(train = cvlist[[l]]$train,
              test = cvlist[[l]]$test, method = "cluster", cvoption = cvopt,
              ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, ...)
          list("minfunval" = funval, "bestnclust" = NA_integer_)
        }
      })
      allObj <- vapply(myObjFun, function(x) x$minfunval, 0)
      bestObj <- min(allObj)
      if (bestObj > bestObj.smallermodel) {
        break # best model of this size is worse than best model of smaller size
      } else {
        bestObj.smallermodel <- bestObj
        bestcovariates <- which(thissubsets[which.min(allObj), ])
        bestnclust <- myObjFun[[which.min(allObj)]]$bestnclust
      }
    }
    result <- list("selectedcols" = bestcovariates,
                   "bestnclust" = bestnclust)
  } else if (!greedy) { # exhaustive search for best subset

    cvlist <- lapply(1:nrow(allsubsets), function(l) {
      fulldat$Xreduced <- fulldat$X[, allsubsets[l, ], drop = FALSE]
      train <- traincalc(fulldat = fulldat, ind = cvind$train,
        method = "cluster", pentype = "L2", nclustmetric = nclust2pass, ...)
      test <- testcalc(fulldat = fulldat, ind = cvind$test,
        method = "cluster", pentype = "L2", cvoption = cvopt,
        train = train, ...)
      list("train" = train, "test" = test)
    })

    myObjFun <- lapply(1:nrow(allsubsets), function(l) {
      if (nclust == "foldcv") {
        if (cvopt == "partitionres")
          stop("This configuration not implemented in varsel.cv.cluster")
        maxnclust <- length(cvlist[[l]]$train[[1]]$allclust)
        allCVObjFunvals <- rep(Inf, maxnclust)
        for (r in 2:maxnclust) {
          allCVObjFunvals[r] <- CVObjFun.varsel(train = cvlist[[l]]$train,
            test = cvlist[[l]]$test, method = "cluster", cvoption = cvopt,
            ctol = 1e-10, solver = "quadprog", Dmatfull = Dm,
            nclust = r, ...)
          if (r >= 4) {
            if (allCVObjFunvals[r] > allCVObjFunvals[r - 1] &&
                allCVObjFunvals[r] > allCVObjFunvals[r - 2] &&
                allCVObjFunvals[r] > allCVObjFunvals[r - 3]) break
          }
        }
        list("minfunval" = min(allCVObjFunvals),
             "bestnclust" = which.min(allCVObjFunvals))
      } else {
        Dm <- if (cvopt == "testsetols") {
          NULL
        } else if (cvopt == "partitionres") {
          clustering <- doclust(X = fulldat$X[, allsubsets[l, ],
            drop = FALSE], nclust = nclust, ...)
          fulldat$L <- matrix(data = 0, nrow = nrow(fulldat$Msq),
            ncol = length(clustering$s))
          for (i in 1:nrow(fulldat$Msq)) {
            fulldat$L[i, clustering$C[i]] <- 1
          }
          fulldat$Msq %*% fulldat$L
        }
        funval <- CVObjFun.varsel(train = cvlist[[l]]$train,
            test = cvlist[[l]]$test, method = "cluster", cvoption = cvopt,
            ctol = 1e-10, solver = "quadprog", Dmatfull = Dm, ...)
        list("minfunval" = funval, "bestnclust" = NA_integer_)
      }
    })
    allObj <- vapply(myObjFun, function(x) x$minfunval, 0)
    bestl <- which.min(allObj)
    bestnclust <- ifelse(nclust == "foldcv",
            myObjFun[[which.min(allObj)]]$bestnclust, NA_integer_)
    result <- list("selectedcols" = which(allsubsets[bestl, ]),
                   "bestnclust" = bestnclust)
    if (nullmodelval <= min(allObj)) { # null model is best
      result <- list("selectedcols" = 1L,
                     "bestnclust" = 1L)
    }
  }
  result
}


# optimise nclust (r) parameter using cross-validation
tune.nclust <- function(fulldat, nfolds = 5L, nclust2search = 1L:20L,
  ctol = 1e-10, solver = "quadprog", cvopt = "testsetols", ...) {

  n <- length(fulldat$esq)

  cvind <- makefolds(fulldat$esq, K = nfolds)
  if (is.null(nclust2search))
    nclust2search <- 1:min(vapply(cvind$train, length, 0L))

  trainvars <- traincalc(fulldat = fulldat, ind = cvind$train,
      method = "cluster", pentype = "L2", ...)
  testvars <- testcalc(fulldat = fulldat, ind = cvind$test,
      method = "cluster", pentype = "L2", cvoption = cvopt,
      train = trainvars, ...)

  myObjFun <- rep(Inf, length(nclust2search))

  if (cvopt == "partitionres") {
    allclust <- lapply(nclust2search, function(r) {
      doclustmult(X = fulldat$Xreduced, nclust = r, ...)$clustering
    })
    allclustbynclust <- NULL
  }

  for (l in seq_along(nclust2search)) {
    if (cvopt == "partitionres") {
      L <- matrix(data = 0, nrow = n, ncol = nclust2search[l])
      for (j in 1:nclust2search[l]) {
        L[1:n %in% allclust[[l]]$clustering$s[[j]], j] <- 1
      }
      passDmat <- fulldat$Msq %*% L
    } else {
      passDmat <- NULL
    }
    myObjFun[l] <- CVObjFun.nclust(nclus = nclust2search[l],
        train = trainvars, test = testvars, solver = solver,
        cvoption = cvopt, Dmatfull = passDmat, ...)

    if (l >= 3) {
      if (myObjFun[l] > myObjFun[l - 1] &
          ((myObjFun[l] - min(myObjFun[1:l])) / min(myObjFun[1:l]) > 1.0)) {
        break # Truncate search as MSEs are increasing
        # and there is a 100% increase compared to best MSE
      }
    }
  }
  nclust2search[which.min(myObjFun)]
}

# mats contains D, Dcross, cvec, A, P, bvec
tune.lambda.qgcv <- function(mats, solver = "osqp",
                             unconstrained = FALSE, ...) {

  # solver must be either osqp or quadprog
  qgcvloss <- function(lam, m = mats, ...) {
    n <- length(m$esq)
    qpsolu <- qpest(Dcross = m$Dcross, P = m$P, A = m$A,
                    cvec = m$cvec, bvec = m$bvec, lambda = lam,
                    solver = solver, returndual = TRUE)
    R <- which(qpsolu$dual != 0)
    Hlam <- m$D %*% solve(m$Dcross + lam * m$P) %*% t(m$D)
    if (length(R) == 0 || unconstrained) {
      tr <- sum(diag(Hlam))
    } else {
      AR <- m$A[R, , drop = FALSE]
      Ulam <- m$D %*% solve(m$Dcross + lam * m$P) %*% t(AR)
      Glam <- handle.error(Ulam %*% solve(AR %*%
        solve(m$Dcross + lam * m$P) %*% t(AR)) %*% t(Ulam))
      tr <- sum(diag(Hlam - Glam))
    }
    esqpred <- drop(m$D %*% qpsolu$coef.est)
    if (length(tr) == 0 || is.null(tr) || is.nan(tr) || is.na(tr)) {
      Inf
    } else {
      sum(((m$esq - esqpred) / (1 - tr / n)) ^ 2) / n
    }
  }

  lamexpseq <- c(10 ^ seq(-5, -1, 1), seq(-0.5, 5, 0.25))
  qgcvlossvals <- vapply(lamexpseq, qgcvloss, 0)
  minind <- which.min(qgcvlossvals)[1]
  if (minind == 1) {
    GSS(f = qgcvloss, a = 0, b = lamexpseq[1], ...)$argmin
  } else if (minind == length(lamexpseq)) {
    GSS(f = qgcvloss, a = lamexpseq[length(lamexpseq)],
        b = 10 ^ 6, ...)$argmin
  } else {
    GSS(f = qgcvloss, a = lamexpseq[minind - 1],
        b = lamexpseq[minind + 1], ...)$argmin
  }
}

# mats contains D, Dcross, cvec, A, P, bvec
# P must be a vector, not matrix, and solver must be roi
tune.lambda.qgcv.lasso <- function(mats, unconstrained = FALSE,
        tr_approx = FALSE, Rtol = 1e-8, ...) {

  qgcvloss <- function(lam, m = mats, appr = tr_approx) {
    n <- length(m$esq)
    # Fit unconstrained model
    qpsol.con <- qpest(Dcross = m$Dcross, P = m$P, A = m$A,
                       cvec = m$cvec, bvec = m$bvec,
                       lambda = lam, solver = "roi")

    if (appr) {
      tr <- sum(abs(qpsol.con$coef.est[-1]) > Rtol)
    } else {
      qpsol.un <- qpest(Dcross = m$Dcross, P = m$P, A = m$A,
        cvec = m$cvec, bvec = rep(-1e6, length(m$bvec)),
        lambda = lam, solver = "roi")

      W <- diag(c(0, abs(qpsol.un$coef.est[-1])))
      Wminus <- MASS::ginv(W)
      R <- which(abs(drop(m$A %*% qpsol.con$coef.est - m$bvec)) <= Rtol)

      Hlam <- m$D %*% solve(m$Dcross + lam * Wminus) %*% t(m$D)
      if (length(R) == 0 || unconstrained) {
        tr <- sum(diag(Hlam))
      } else {
        AR <- m$A[R, , drop = FALSE]
        Ulam <- m$D %*% solve(m$Dcross + lam * Wminus) %*% t(AR)
        Glam <- handle.error(Ulam %*% Rfast::spdinv(AR %*%
          solve(m$Dcross + lam * Wminus) %*% t(AR)) %*% t(Ulam))
        tr <- sum(diag(Hlam - Glam))
      }
    }
    esqpred <- drop(m$D %*% qpsol.con$coef.est)
    if (length(tr) == 0 || is.null(tr) || is.nan(tr) || is.na(tr)) {
      Inf
    } else {
      sum(((m$esq - esqpred) / (1 - tr / n)) ^ 2) / n
    }
  }

  lamexpseq <- c(10 ^ seq(-5, -1, 1), seq(-0.5, 5, 0.25))
  qgcvlossvals <- vapply(lamexpseq, qgcvloss, 0)
  minind <- which.min(qgcvlossvals)[1]
  if (minind == 1) {
    GSS(f = qgcvloss, a = 0, b = lamexpseq[1], ...)$argmin
  } else if (minind == length(lamexpseq)) {
    GSS(f = qgcvloss, a = lamexpseq[length(lamexpseq)],
        b = 10 ^ 6, ...)$argmin
  } else {
    GSS(f = qgcvloss, a = lamexpseq[minind - 1],
        b = lamexpseq[minind + 1], ...)$argmin
  }
}

tune.lambda.cv <- function(fulldat, method = NULL,
      nfolds = 5L, d = 2L, lambdaUB = 1e5, solver = "quadprog",
      pentype = "L2", cvopt = "testsetols", tsk = NULL, tsm = NULL,
      Dmat = NULL, Pmat = NULL, ...) {

    cvind <- makefolds(fulldat$esq, K = nfolds)
    trainvars <- traincalc(fulldat = fulldat, ind = cvind$train,
        method = method, pentype = pentype,
        tsk = tsk, tsm = tsm, polydeg = d, ...)

    testvars <- testcalc(fulldat = fulldat, ind = cvind$test,
        method = method, pentype = pentype, cvoption = cvopt,
        polydeg = d, ...)

    myf <- function(x, ...) CVObjFun.lambda(lam = x, train = trainvars,
      test = testvars, meth = method, solver = solver, cvoption = cvopt,
      Dmatfull = Dmat, ...)

  thex0ub <- lambdaUB
  if (method == "polynomial" && pentype == "L1") {
    # Find smallest upper bound where all coefficient estimates
    #  are 0 other than intercept
    for (mm in log10(lambdaUB):-3) {
      thiscoef <- getL1coef(lam = 10 ^ mm, train = trainvars,
          test = testvars, solver = solver, ...)
      if (any(abs(thiscoef[-1]) > 1e-10)) {
        break
      }
    }
    thex0ub <- 10 ^ mm
  }

  ptstocheck <- 10 ^ (-3:log10(thex0ub))
  trainerrscheck <- vapply(ptstocheck, myf, 0)

  finalub <- ptstocheck[which.min(trainerrscheck) + 1]
  if (length(finalub) == 0L) stop("!!!!!")
  if (is.na(finalub)) finalub <- ptstocheck[which.min(trainerrscheck)]
  br <- bracket(f = myf, x0 = finalub, ...)
  mygss <- GSS(f = myf, a = br$L, b = br$U, ...)
  mygss$argmin
}

CVObjFun.varsel <- function(train, test, method, cvoption = "testsetols",
                            ctol = 1e-10, solver = "quadprog",
                            Dmatfull = NULL, nclust = NULL, ...) {

  msqErr <- mean(vapply(1:length(train), function(j) {
    if (method == "linear") {
      mats <- train[[j]]
    } else if (method == "cluster") {
      if (is.null(nclust)) {
        mats <- train[[j]]$onemats
      } else {
        mats <- train[[j]]$allmats[[nclust]]
      }
    }
    qpsol <- qpest(Dcross = mats$Dcross, P = mats$P, A = mats$A,
      cvec = mats$cvec, bvec = rep(ctol, nrow(mats$A)), lambda = 0,
      solver = solver, ...)
    if (cvoption == "testsetols") {
      if (method == "linear") {
        deviation <- test[[j]]$esq - drop(test[[j]]$D %*% qpsol$coef.est)
      } else if (method == "cluster") {
        if (is.null(nclust)) {
          deviation <- test[[j]]$esq - drop(test[[j]]$Dmat %*% qpsol$coef.est)
        } else {
          deviation <- test[[j]]$esq -
            drop(test[[j]]$allDmats[[nclust]] %*% qpsol$coef.est)
        }
      }
    } else if (cvoption == "partitionres") {
      allesqpred <- drop(Dmatfull %*% qpsol$coef.est)
      deviation <- test[[j]]$esq - allesqpred[test[[j]]$ind]
    }
    drop(crossprod(deviation)) / length(test[[j]]$esq)
  }, 0))
  msqErr
}

# Finds cross-validated mean of objective function
#   for given value of nclust clustering model
CVObjFun.nclust <- function(nclus, train, test, cvoption = "testsetols",
      ctol = 1e-10, solver = "quadprog",
      Dmatfull = NULL, ...) {

  if ("omega.hat" %in% names(train[[1]])) {
    msqErr <- mean(vapply(1:length(train), function(j) {
      sum((test[[j]]$esq - train[[j]]$omega.hat) ^ 2) / length(test[[j]]$esq)
    }, 0))
  } else {
    msqErr <- mean(vapply(1:length(train), function(j) {
      qpsol <- qpest(Dcross = train[[j]]$allmats[[nclus]]$Dcross,
                     P = train[[j]]$allmats[[nclus]]$P,
                     A = train[[j]]$allmats[[nclus]]$A,
                     cvec = train[[j]]$allmats[[nclus]]$cvec,
                     bvec = rep(ctol, nclus), solver = solver, ...)
      if (cvoption == "testsetols") {
        deviation <- test[[j]]$esq -
          drop(test[[j]]$allDmats[[nclus]] %*% qpsol$coef.est)
      } else if (cvoption == "partitionres") {
        allesqpred <- drop(Dmatfull %*% qpsol$coef.est)
        deviation <- test[[j]]$esq - allesqpred[test[[j]]$ind]
      }
      drop(crossprod(deviation)) / length(test[[j]]$esq)
    }, 0))
  }
  msqErr
}

getL1coef <- function(lam, train, test, solver = "quadprog",
                      ctol = 1e-10, ...) {
  j <- 1
  qpsol <- qpest(Dcross = train[[j]]$Dcross, P = train[[j]]$P,
                 A = train[[j]]$A, cvec = train[[j]]$cvec,
                 bvec = rep(ctol, nrow(train[[j]]$A)), lambda = lam,
                 solver = solver, ...)
  qpsol$coef.est
}

# Finds cross-validated MSE for a given lambda
CVObjFun.lambda <- function(lam, train, test, meth, ctol = 1e-10,
    solver = "quadprog", cvoption = "testsetols", Dmatfull = NULL, ...) {

  if ("omega.hat" %in% names(train[[1]])) {
    msqErr <- mean(vapply(1:length(train), function(j) {
      sum((test[[j]]$esq - train[[j]]$omega.hat) ^ 2) /
        length(test[[j]]$esq)
    }, 0))
  } else {
    msqErr <- mean(vapply(1:length(train), function(j) {
      qpsol <- qpest(Dcross = train[[j]]$Dcross, P = train[[j]]$P,
        A = train[[j]]$A, cvec = train[[j]]$cvec,
        bvec = rep(ctol, nrow(train[[j]]$A)), lambda = lam,
        solver = solver, ...)
      if (cvoption == "testsetols") {
        if (meth == "polynomial") {
          deviation <- test[[j]]$esq - drop(test[[j]]$D %*% qpsol$coef.est)
        } else if (meth == "spline") {
          train[[j]]$gam.thin$coefficients <- qpsol$coef.est
          newdf <- data.frame("esq" = test[[j]]$esq, "X" = test[[j]]$Xreduced)
          Dtest <- cbind(1, mgcv::PredictMat(train[[j]]$gam.thin$smooth[[1]],
                      data = newdf))
          deviation <- test[[j]]$esq - drop(Dtest %*% qpsol$coef.est)
        }
      } else if (cvoption == "partitionres") {
        allesqpred <- drop(Dmatfull %*% qpsol$coef.est)
        deviation <- test[[j]]$esq - allesqpred[test[[j]]$ind]
      }
      drop(crossprod(deviation)) / length(test[[j]]$esq)
    }, 0))
  }
  msqErr
}

# for cluster model with only 1 cluster
traincalc0 <- function(fulldat, ind) {
  lapply(1:length(ind), function(j) {
    Xfold <- fulldat$X[ind[[j]], , drop = FALSE]
    yfold <- fulldat$y[ind[[j]]]
    Msqfold <- (diag(length(yfold)) -
                  Xfold %*% Rfast::spdinv(crossprod(Xfold)) %*% t(Xfold)) ^ 2
    esqfold <- stats::resid(stats::lm.fit(x = Xfold, y = yfold)) ^ 2
    Lfold <- matrix(data = 1L, nrow = nrow(Xfold), ncol = 1L)
    Dfold <- Msqfold %*% Lfold
    Afold <- diag(1)
    Pfold <- matrix(data = 0, nrow = 1L, ncol = 1L)

    matsfold <- list("D" = Dfold, "Dcross" = crossprod(Dfold),
                     "cvec" = t(Dfold) %*% esqfold,
                     "A" = Afold, "P" = Pfold)
    return(list("allclust" = NA, "allmats" = list(matsfold),
                "Xreduced" = NA, "esq" = esqfold))
  })
}

testcalc0 <- function(fulldat, ind, cvoption = "testsetols") {
  lapply(1:length(ind), function(j) {
    Xfold <- fulldat$X[ind[[j]], , drop = FALSE]
    yfold <- fulldat$y[ind[[j]]]
    if (cvoption == "testsetols") {
      Msqfold <- (diag(length(yfold)) - Xfold %*%
                    Rfast::spdinv(crossprod(Xfold)) %*% t(Xfold)) ^ 2
      esqfold <- stats::resid(stats::lm.fit(x = Xfold, y = yfold)) ^ 2
      Lfold <- matrix(data = 1, nrow = nrow(Xfold), ncol = 1)
      Dfold <- Msqfold %*% Lfold
      return(list("allDmats" = list(Dfold), "esq" = esqfold))
    } else if (cvoption == "partitionres") {
      return(list("esq" = fulldat$esq[ind[[j]]], "ind" = ind[[j]]))
    }
  })
}

# Performs calculations on training sets
traincalc <- function(fulldat, ind, method, pentype = "L2",
                      tsk = NULL, tsm = NULL,
                      nclustmetric = NULL, maxnclust = 20L,
                      polydeg = 2L, ...) {

  lapply(1:length(ind), function(j) {
    Xfold <- fulldat$X[ind[[j]], , drop = FALSE]
    yfold <- fulldat$y[ind[[j]]]
    Msqfold <- (diag(length(yfold)) -
        Xfold %*% Rfast::spdinv(crossprod(Xfold)) %*% t(Xfold)) ^ 2
    esqfold <- stats::resid(stats::lm.fit(x = Xfold, y = yfold)) ^ 2
    Xreducedfold <- fulldat$Xreduced[ind[[j]], , drop = FALSE]
    colnames(Xreducedfold) <- NULL

    if (method == "linear") {
      Lfold <- Afold <- Xreducedfold
      Pfold <- matrix(data = 0, nrow = ncol(Lfold), ncol = ncol(Lfold))
    } else if (method == "polynomial") {
      Lfold <- Afold <- makepolydesign(Xreducedfold, polydeg)
      if (pentype == "L2") {
        Pfold <- diag(ncol(Lfold))
        Pfold[1, 1] <- 0
      } else if (pentype == "L1") {
        Pfold <- rep(1, 2 * ncol(Lfold))
        Pfold[c(1, ncol(Lfold) + 1)] <- 0
      }
    } else if (method == "spline") {
      thindat <- data.frame("esq" = esqfold, "X" = Xreducedfold)
      if (!is.null(tsk) && !is.null(tsm)) {
        argbits <- paste0(", k = ", tsk, ", m = ", tsm)
      } else if (!is.null(tsk) && is.null(tsm)) {
        argbits <- paste0(", k = ", tsk)
      } else if (is.null(tsk) && !is.null(tsm)) {
        argbits <- paste0(", m = ", tsm)
      } else {
        argbits <- character(0)
      }
      if (ncol(Xreducedfold) > 1) {
        thinformula <- stats::as.formula(paste0("esq ~ s(",
          paste0("X.", 1:ncol(Xreducedfold), collapse = ","),
          ", bs = \"tp\"", argbits, ")"))
      } else {
        thinformula <- stats::as.formula(paste0("esq ~ s(X, bs = \"tp\"",
                                         argbits, ")"))
      }
      gam.thin <- mgcv::gam(thinformula, data = thindat, ...)
      Dfold <- cbind(1, mgcv::PredictMat(gam.thin$smooth[[1]],
                                         data = thindat))
      Lfold <- Afold <- Rfast::spdinv(Msqfold) %*% Dfold

      sobj <- eval(parse(text = as.character(thinformula)[3]))
      Pfold <- mgcv::smooth.construct.tp.smooth.spec(sobj,
                data = thindat, knots = NULL)$S[[1]]
      return(list("D" = Dfold, "Dcross" = crossprod(Dfold),
                  "cvec" = t(Dfold) %*% esqfold,
                  "A" = Afold, "P" = Pfold, "gam.thin" = gam.thin))
    } else if (method == "cluster") {
      if (is.null(nclustmetric) || nclustmetric == "foldcv") {
        nclustupperlim <- min(nrow(Xreducedfold), maxnclust)
        allclust <- doclustmult(X = Xreducedfold,
                                nclust = 1L:nclustupperlim, ...)
        matsfold <- lapply(1:nclustupperlim, function(r) {
          L <- matrix(data = 0, nrow = nrow(Xreducedfold),
                      ncol = allclust$clustering[[r]]$nclust)
          for (i in 1:nrow(Xreducedfold)) {
            L[i, allclust$clustering[[r]]$C[i]] <- 1
          }
          D <- Msqfold %*% L
          A <- diag(allclust$clustering[[r]]$nclust)
          P <- matrix(data = 0, nrow = allclust$clustering[[r]]$nclust,
                      ncol = allclust$clustering[[r]]$nclust)
          list("D" = D, "Dcross" = crossprod(D),
               "cvec" = t(D) %*% esqfold,
               "A" = A, "P" = P)
        })
        return(list("allclust" = allclust, "allmats" = matsfold,
                    "Xreduced" = Xreducedfold, "esq" = esqfold))
      } else {
        oneclust <- doclust(X = Xreducedfold, nclust = nclustmetric, ...)
        Lfold <- matrix(data = 0, nrow = nrow(Xreducedfold),
                        ncol = oneclust$nclust)
        for (i in 1:nrow(Xreducedfold)) {
          Lfold[i, oneclust$C[i]] <- 1
        }
        Dfold <- Msqfold %*% Lfold
        Afold <- diag(oneclust$nclust)
        Pfold <- matrix(data = 0, nrow = oneclust$nclust,
                        ncol = oneclust$nclust)
        onemats <- list("D" = Dfold, "Dcross" = crossprod(Dfold),
                        "cvec" = t(Dfold) %*% esqfold,
                        "A" = Afold, "P" = Pfold)
        return(list("oneclust" = oneclust, "onemats" = onemats,
                    "Xreduced" = Xreducedfold, "esq" = esqfold))
      }
    }
    Dfold <- Msqfold %*% Lfold
    return(list("Dcross" = crossprod(Dfold),
                "cvec" = t(Dfold) %*% esqfold, "A" = Afold,
                "P" = Pfold))
  })
}

testcalc <- function(fulldat, ind, method, pentype = "L2",
          cvoption = "testsetols", train = NULL, polydeg = 2L, ...) {

  lapply(1:length(ind), function(j) {
    Xreducedfold <- fulldat$Xreduced[ind[[j]], , drop = FALSE]
    Xfold <- fulldat$X[ind[[j]], , drop = FALSE]
    yfold <- fulldat$y[ind[[j]]]
    if (cvoption == "testsetols") {
      Msqfold <- (diag(length(yfold)) - Xfold %*%
            Rfast::spdinv(crossprod(Xfold)) %*% t(Xfold)) ^ 2
      esqfold <- stats::resid(stats::lm.fit(x = Xfold, y = yfold)) ^ 2

      if (method == "linear") {
        Lfold <- Xreducedfold
        Pfold <- matrix(data = 0, nrow = ncol(Lfold), ncol = ncol(Lfold))
      } else if (method == "polynomial") {
        Lfold <- makepolydesign(Xreducedfold, polydeg)
        if (pentype == "L2") {
          Pfold <- diag(ncol(Lfold))
          Pfold[1, 1] <- 0
        } else if (pentype == "L1") {
          Pfold <- rep(1, 2 * ncol(Lfold))
          Pfold[c(1, ncol(Lfold) + 1)] <- 0
        }
      } else if (method == "cluster") {
        if ("allclust" %in% names(train[[j]])) {
          # clustering was done for all numbers of clusters from 1:maxnclust
          allDmats <- lapply(1:length(train[[j]]$allclust$clustering),
            function(r) {
              testC <- add2clust(old = train[[j]]$Xreduced, new = Xreducedfold,
                s = train[[j]]$allclust$clustering[[r]]$s,
                distmetric = train[[j]]$allclust$distmetric,
                linkage = train[[j]]$allclust$linkage,
                scaledata = train[[j]]$allclust$scaledata)
              testL <- matrix(data = 0, nrow = nrow(Xreducedfold),
                ncol = train[[j]]$allclust$clustering[[r]]$nclust)
              for (i in 1:nrow(Xreducedfold)) {
                testL[i, testC[i]] <- 1
              }
              Msqfold %*% testL
            })
          return(list("allDmats" = allDmats, "Xreduced" = Xreducedfold,
                      "esq" = esqfold))
        } else if ("oneclust" %in% names(train[[j]])) {
          # clustering was done for only one no. of clusters, from elbow method
          testC <- add2clust(old = train[[j]]$Xreduced, new = Xreducedfold,
                             s = train[[j]]$oneclust$s,
                             distmetric = train[[j]]$oneclust$distmetric,
                             linkage = train[[j]]$oneclust$linkage,
                             scaledata = train[[j]]$oneclust$scaledata)
          testL <- matrix(data = 0, nrow = nrow(Xreducedfold),
                          ncol = train[[j]]$oneclust$nclust)
          for (i in 1:nrow(Xreducedfold)) {
            testL[i, testC[i]] <- 1
          }
          Dmat <- Msqfold %*% testL
          return(list("Dmat" = Dmat, "Xreduced" = Xreducedfold,
                      "esq" = esqfold))
        }
      } else if (method == "spline") {
        # Can't compute other matrices in these cases.
        return(list("esq" = esqfold, "Msq" = Msqfold,
                    "Xreduced" = Xreducedfold))
      }
      Dfold <- Msqfold %*% Lfold
      return(list("D" = Dfold, "Xreduced" = Xreducedfold, "esq" = esqfold))
    } else if (cvoption == "partitionres") {
      return(list("esq" = fulldat$esq[ind[[j]]], "ind" = ind[[j]],
                  "Xreduced" = Xreducedfold))
    }
  })
}

makefolds <- function(esq, K = 5L) {
  testfolds <- caret::createFolds(esq, K)
  trainfolds <- lapply(testfolds, function(x) setdiff(1:length(esq), x))
  list("train" = trainfolds, "test" = testfolds)
}

# P will be either a q x q identity matrix (L2) with [1, 1] element = 0
#  or a 2q-vector of ones (L1) with 1st and (q+1)st elements = 0
qpest <- function(Dcross, P, A, cvec, bvec, lambda = 0,
                  solver = c("quadprog", "quadprogXT",
                             "roi", "osqp"),
                  returndual = FALSE, ...) {

  q <- length(cvec)
  coef.hat.pm <- NULL
  dual <- NULL
  pentype <- ifelse(is.matrix(P), "L2", "L1")
  solver <- match.arg(solver, c("quadprog", "quadprogXT",
                                "roi", "osqp"))
  if (solver == "quadprog") {
    qpfunc <- quadprog::solve.QP.compact
    mycompact <- TRUE
    if (pentype == "L2") {
      Qmat <- Dcross + lambda * P
      Qeigen <- eigen(Qmat, only.values = TRUE)$values
      if (any(Qeigen < 1e-10)) {
        Qmat <- as.matrix(Matrix::nearPD(Qmat,
                  eig.tol = 1e-10, posd.tol = 1e-10)$mat)
      }
      qp <- quadprogXT::buildQP(Dmat = Qmat,
              dvec = cvec, Amat = t(A),
              bvec = bvec, dvecPosNeg = NULL, compact = mycompact)
      qpsolu <- do.call(what = qpfunc, args = qp)
      coef.hat <- qpsolu$solution
      if (returndual) dual <- qpsolu$Lagrangian
    } else if (pentype == "L1") {
      newDcross <- rbind(cbind(Dcross, -Dcross), cbind(-Dcross, Dcross))
      newDcross2 <- as.matrix(Matrix::nearPD(newDcross,
              eig.tol = 1e-10, posd.tol = 1e-10)$mat)
      newcvec <- c(cvec, -cvec)
      newA <- rbind(cbind(A, -A), diag(2 * q))
      newbvec <- c(bvec, rep(0, 2 * q))
      qp <- quadprogXT::buildQP(Dmat = newDcross2,
              dvec = newcvec - lambda / 2 * P, Amat = t(newA),
              bvec = newbvec, dvecPosNeg = NULL, compact = mycompact)
      qpsolu <- do.call(what = qpfunc, args = qp)
      coef.hat.pm <- qpsolu$solution
      coef.hat <- coef.hat.pm[1:q] - coef.hat.pm[(q + 1):(2 * q)]
      if (returndual) dual <- qpsolu$Lagrangian
    }
  } else if (solver == "quadprogXT") {
    qpfunc <- quadprog::solve.QP
    if (pentype == "L2") {
      qp <- quadprogXT::buildQP(Dmat = Dcross + lambda * P,
              dvec = cvec, Amat = t(A), bvec = bvec,
              dvecPosNeg = NULL, compact = TRUE)
    } else if (pentype == "L1") {
      qp <- quadprogXT::buildQP(Dmat = Dcross,
              dvec = cvec, Amat = t(A), bvec = bvec,
              dvecPosNeg = -lambda / 2 * P, compact = TRUE)
    }
    coef.hat.pm <- do.call(what = qpfunc, args = qp)$solution
    coef.hat <- coef.hat.pm[1:q]
  } else if (solver == "roi") {
    if (pentype == "L2") {
      Qsparse <- slam::as.simple_triplet_matrix(Dcross + lambda * P)
      Asparse <- slam::as.simple_triplet_matrix(A)
      qpsparse <- ROI::OP(ROI::Q_objective(Q = Qsparse, L = -cvec),
            ROI::L_constraint(L = Asparse,
            dir = rep(">=", length(bvec)), rhs = bvec))
      coef.hat <- ROI::ROI_solve(qpsparse, solver = "qpoases")$solution
    } else if (pentype == "L1") {
      newDcross <- rbind(cbind(Dcross, -Dcross), cbind(-Dcross, Dcross))
      newcvec <- c(cvec, -cvec)
      newA <- rbind(cbind(A, -A), diag(2 * q))
      newbvec <- c(bvec, rep(0, 2 * q))
      Qsparse <- slam::as.simple_triplet_matrix(newDcross)
      Asparse <- slam::as.simple_triplet_matrix(newA)
      qpsparse <- ROI::OP(ROI::Q_objective(Q = Qsparse,
            L = -newcvec + lambda / 2 * P),
            ROI::L_constraint(L = Asparse,
            dir = rep(">=", length(newbvec)), rhs = newbvec))
      coef.hat.pm <- ROI::ROI_solve(qpsparse, solver = "qpoases")$solution
      coef.hat <- coef.hat.pm[1:q] - coef.hat.pm[(q + 1):(2 * q)]
    }
  } else if (solver == "osqp") {
    if (pentype == "L2") {
      Qmat <- Dcross + lambda * P
      Qeigen <- eigen(Qmat, only.values = TRUE)$values
      if (any(Qeigen < 1e-10)) {
        Qmat <- as.matrix(Matrix::nearPD(Qmat,
              eig.tol = 1e-10, posd.tol = 1e-10)$mat)
      }
      qpsolu <- osqp::solve_osqp(P = Qmat,
              q = -cvec, A = A, l = bvec,
              pars = osqp::osqpSettings(verbose = FALSE, polish = TRUE))
      coef.hat <- qpsolu$x
      if (returndual) dual <- -qpsolu$y
    } else if (pentype == "L1") {
      newDcross <- rbind(cbind(Dcross, -Dcross), cbind(-Dcross, Dcross))
      newcvec <- c(cvec, -cvec)
      newA <- rbind(cbind(A, -A), diag(2 * q))
      newbvec <- c(bvec, rep(0, 2 * q))
      qpsolu <- osqp::solve_osqp(P = newDcross,
          q = -newcvec + lambda / 2 * P, A = newA, l = newbvec,
          pars = osqp::osqpSettings(verbose = FALSE, polish = TRUE))
      coef.hat.pm <- qpsolu$x
      coef.hat <- coef.hat.pm[1:q] - coef.hat.pm[(q + 1):(2 * q)]
      if (returndual) dual <- -qpsolu$y
    }
  }
  list("coef.est" = coef.hat, "coef.est.pm" = coef.hat.pm,
       "dual" = dual)
}

makepolydesign <- function(X, d = 2L) {
  dimm <- dim(X)
  n <- dimm[1]
  p <- dimm[2]

  # Each row of polyind2 contains column indices of X in one term
  # of d-degree polynomial
  # where 1's represent the intercept column
  polyind <- as.matrix(expand.grid(lapply(1:d, function(j) 1:p)))
  polyind2 <- unique(Rfast::rowSort(polyind))

  thispolyX <- as.matrix(as.data.frame(lapply(1:nrow(polyind2),
    function(m) {
      vapply(1:n, function(i) prod(X[i, polyind2[m, ]]), 0)
    })))

  colnames(thispolyX) <- paste0("V", 1:ncol(thispolyX))
  thispolyX[, !duplicated(t(thispolyX)), drop = FALSE]
}

# Uses bracketing method (See Wu & Lange 2008) to narrow down
#  interval for golden section search
bracket <- function(f, x0 = 1e5, rr = 0.5, maxbrack = 25L,
                    minchange = 1e-4, ...) {

  if (rr <= 0 || rr >= 1) stop("rr must be in interval (0,1)")
  xseq <- x0 * rr ^ (0:maxbrack)
  prevfval <- f(xseq[1], ...)
  thisfval <- f(xseq[2], ...)
  foundit <- FALSE
  for (j in 3:(maxbrack + 1)) {
    nextfval <- f(xseq[j], ...)
    if ((thisfval - prevfval < -minchange) &&
        (thisfval - nextfval < -minchange)) {
      foundit <- TRUE
      break
    }
    prevfval <- thisfval
    thisfval <- nextfval
  }
  if (foundit) {
    list("U" = xseq[j - 2], "L" = xseq[j], "iterations" = j)
  } else {
    warning("Bracketing failed to find appropriate interval. Upper bound for golden section search will be set to initial upper bound.")
    list("U" = x0, "L" = 0, "iterations" = j)
  }
}

doclust <- function(X, distmetric = "euclidean",
                    linkage = "complete", scaledata = TRUE,
                    nclust = c("elbow.swd", "elbow.mwd",
                               "elbow.both"),
                    returnmetric = FALSE, ...) {

  hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) X <- X[, -hasintercept[[2]], drop = FALSE]
  if (scaledata) X <- scale(X)
  n <- nrow(X)
  metric <- NA

  dm <- stats::dist(X, method = distmetric)
  agglom <- stats::hclust(dm, method = linkage)

  if (is.numeric(nclust)) {
    myr <- nclust
    if (returnmetric) metric <- list()
  } else {
    nclust <- match.arg(nclust, c("elbow.swd", "elbow.mwd",
                                  "elbow.both"))
    if (nclust == "elbow.both") {

      doswd <- SWDelbow(agglom, X = X, dmet = distmetric)
      domwd <- MWDelbow(agglom)
      myr <- round(mean(c(doswd$r, domwd$r)))

      if (returnmetric) metric <- list("mwd" = domwd$mwd, "swd" = doswd$swd)

    } else if (nclust == "elbow.mwd") {
      domwd <- MWDelbow(agglom)
      myr <- domwd$r
      if (returnmetric) metric <- list("mwd" = domwd$mwd)
    } else if (nclust == "elbow.swd") {

      doswd <- SWDelbow(agglom, X = X, dmet = distmetric)
      myr <- doswd$r
      if (returnmetric) metric <- list("swd" = doswd$swd)
    }
  }

  clusterind <- stats::cutree(agglom, k = myr)
  clusterlist <- lapply(1:myr, function(i) {
    which(clusterind == i)
  })

  val <- list("s" = clusterlist, "C" = clusterind, "metric" = metric,
              "distmetric" = distmetric, "linkage" = linkage,
              "scaledata" = scaledata, "nclust" = myr)
  class(val) <- "doclust"
  val
}

# Runs doclust for all permutations of selected variables
doclustperm <- function(X, ...) {
  hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) {
    cols <- 2:ncol(X)
  } else {
    cols <- 1:ncol(X)
  }

  mygrid <- as.matrix(expand.grid(lapply(1:length(cols),
              function(x) c(FALSE, TRUE))))
  if (hasintercept[[1]]) {
    mygrid <- cbind(FALSE, mygrid)
  }
  mygrid <- mygrid[-1, , drop = FALSE]

  doclusts <- lapply(1:nrow(mygrid), function(i) {
    doclust(X[, which(mygrid[i, ]), drop = FALSE], ...)
  })

  usedcols <- vapply(1:nrow(mygrid), function(i) {
    paste(which(mygrid[i, ]), collapse = "")
  }, "")
  names(doclusts) <- usedcols

  class(doclusts) <- "doclustperm"
  doclusts
}

# Assigns observations at random to a fixed number of clusters with
# equal probability. The purpose is so that, when computing
# bootstrap CIs, the size of parameter vector is always the same
# for cluster model even if no variables were selected (homoskedas)
randclust <- function(n, nclust) {
  clusterind <- sample(1:nclust, size = n, replace = TRUE)
  clusterlist <- lapply(1:nclust, function(i) {
    which(clusterind == i)
  })
  list("s" = clusterlist, "C" = clusterind, "metric" = list(),
       "distmetric" = "none", "linkage" = "none",
       scaledata = FALSE)
}


# Returns cluster info for each of a range of nclust values
doclustmult <- function(X, nclust = 1L:n, distmetric = "euclidean",
                        linkage = "complete", scaledata = TRUE, ...) {

  if (!is.numeric(nclust)) stop("nclust must be an integer vector")
  hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) X <- X[, -hasintercept[[2]], drop = FALSE]
  if (scaledata) X <- scale(X)
  n <- nrow(X)

  dm <- stats::dist(X, method = distmetric)
  agglom <- stats::hclust(dm, method = linkage)

  results <- lapply(nclust, function(r) {
    clusterind <- stats::cutree(agglom, k = r)
    clusterlist <- lapply(1L:r, function(i) {
      which(clusterind == i)
    })
    list("s" = clusterlist, "C" = clusterind, "nclust" = r)
  })

  val <- list("clustering" = results, "metric" = NA_real_,
              "distmetric" = distmetric, "linkage" = linkage,
              "scaledata" = scaledata)
  class(val) <- "doclustmult"
  val
}

SWDelbow <- function(anhclust, X = X, dmet = "euclidean", ...) {

  swd <- vapply(1:length(anhclust$order), function(i) {
    C <- stats::cutree(anhclust, k = i)
    sum(vapply(1:i, function(j) {
      sthis <- which(C == j)
      if (length(sthis) == 1) {
        0
      } else {
        wd <- as.matrix(stats::dist(X[sthis, , drop = FALSE], method = dmet))
        sum(wd[upper.tri(wd)])
      }
    }, 0))
  }, 0)
  myr <- inflection::uik(1:length(anhclust$order), swd)
  list("r" = myr, "swd" = swd[myr])
}

MWDelbow <- function(anhclust) {
  mwd <- rev(c(0, anhclust$height))
  myr <- inflection::uik(1:length(anhclust$order), mwd)
  list("r" = myr, "mwd" = mwd[myr])
}

# Assign observation(s) to the nearest existing cluster(s)
add2clust <- function(old, new, s = NULL,
                      distmetric = "euclidean",
                      linkage = "complete",
                      scaledata = TRUE, ...) {

  if (is.null(s)) {
    s <- doclust(X = old)$s # uses default nclust
  }

  if (scaledata) {
    both <- scale(rbind(old, new))
    old <- both[1:nrow(old), , drop = FALSE]
    new <- both[-c(1:nrow(old)), , drop = FALSE]
  }

  nnew <- nrow(new)

  # Calculate distance of each point to existing clusters
  if (linkage == "complete") {
    d2clust <- sapply(1:nrow(new), function(i) {
      vapply(1:length(s), function(j) {
        max(as.matrix(stats::dist(rbind(new[i, ],
            old[s[[j]], ]), method = distmetric))[1, -1])
      }, 0)
    }, simplify = "matrix")
  } else if (linkage == "average") {
    d2clust <- sapply(1:nrow(new), function(i) {
      vapply(1:length(s), function(j) {
        mean(as.matrix(stats::dist(rbind(new[i, ],
            old[s[[j]], ]), method = distmetric))[1, -1])
      }, 0)
    }, simplify = "array")
  } else {
    stop("Only complete and average linkage rules supported in add2clust")
  }
  if (is.matrix(d2clust)) {
    d2clust <- t(d2clust)
  } else {
    d2clust <- as.matrix(d2clust)
  }
  whichclust <- vapply(1:nnew, function(i) {
    which.min(d2clust[i, ])
  }, 0L)
  whichclust
}

hetvarsel <- function(mainlm, hasintercept = NULL, testname = "evans_king",
              alpha = 0.1, bpsq = TRUE, selectonlyone = FALSE, ...) {

  if (inherits(mainlm, "lm")) {
    X <- stats::model.matrix(mainlm)
    y <- stats::model.response(stats::model.frame(mainlm))
    e <- stats::resid(mainlm)
    n <- nrow(X)
    p <- ncol(X)
  } else { # other type of list
    processmainlm(m = mainlm)
  }
  testname <- match.arg(testname, c("breusch_pagan", "bamset",
              "carapeto_holt", "evans_king", "goldfeld_quandt",
              "harrison_mccabe", "honda", "horn", "szroeter"))

  if (alpha < 0 || alpha > 1)
    stop("Significance level must be between 0 and 1")
  dots <- list(...)
  dots$mainlm <- mainlm
  if ("deflator" %in% names(dots)) dots$deflator <- NULL
  testfunc <- utils::getFromNamespace(x = testname, ns = "skedastic")
  dots[!(names(dots) %in% names(formals(testfunc)))] <- NULL

  if (is.null(hasintercept))
    hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) { # drop intercept so it isn't part of X
    deflators <- setdiff(1:p, hasintercept[[2]])
  } else {
    deflators <- 1:p
  }

  if (testname == "breusch_pagan") {
    pvals <- vapply(deflators, function(d) {
      thisdots <- dots
      thisdots$auxdesign <- X[, d, drop = FALSE]
      if (bpsq) thisdots$auxdesign <- cbind(thisdots$auxdesign,
                                            X[, d] ^ 2)
      do.call(what = testfunc, args = thisdots)$p.value
    }, 0)
  } else {
    pvals <- vapply(deflators, function(d) {
      thisdots <- dots
      thisdots$deflator <- d
      do.call(what = testfunc, args = thisdots)$p.value
    }, 0)
  }
  if (selectonlyone) {
    if (any(pvals < alpha)) {
      list("selectedcols" = deflators[which.min(pvals)],
           "p.values" = pvals)
    } else {
      list("selectedcols" = integer(0), "p.values" = pvals)
    }
  } else {
    list("selectedcols" = deflators[which(pvals < alpha)],
         "p.values" = pvals)
  }
}

handle.error <- function(expr, returniferror = NULL) {
  expr_name <- deparse(substitute(expr))
  test <- try(expr, silent = TRUE)
  if (inherits(test, "try-error")) {
    returniferror
  } else {
    test
  }
}

is.error <- function(expr) {
  expr_name <- deparse(substitute(expr))
  test <- try(expr, silent = TRUE)
  inherits(test, "try-error")
}

bootavm <- function(object, bootobject,
                    retune = FALSE, returnall = FALSE,
                    rm_on_constraint = FALSE, rm_nonconverged = TRUE,
                    Brequired = floor(0.75 * length(bootobject)),
                    ...) {

  if (!(inherits(object, "alvm.fit") ||
        inherits(object, "anlvm.fit")))
    stop("object should be of class alvm.fit or anlvm.fit")
  if (!inherits(bootobject, "bootlm"))
    stop("bootobject should be of class bootlm")

  B <- length(bootobject)
  if (retune) {
    thevarselect <- object$selectinfo$varselect
  } else {
    thevarselect <- object$selectinfo$selectedcols
  }

  if (inherits(object, "alvm.fit")) {
    thelambda <- ifelse(retune, object$hyperpar$lambdamethod,
                        object$hyperpar$lambda)

    if (rm_on_constraint) {
      bootaux <- vector("list", length = Brequired)
      a <- 0
      b_used <- rep(NA_integer_, Brequired)
      bootauxbad <- vector("list", length = B)

      if (object$method == "spline") {
        boot_tsk <- object$hyperpar$tsk / 2
      } else {
        boot_tsk <- NULL
      }
      for (b in 1:B) {
        thisfit <- alvm.fit(mainlm = bootobject[[b]],
            varselect = thevarselect,
            model = object$method, lambda = thelambda,
            nclust = object$hyperpar$nclustmethod,
            polypen = object$pentype, solver = object$solver,
            tsk = boot_tsk, reduce2homosked = FALSE, ...)
        if (all(thisfit$var.est > 1.1 * object$constol)) {
          a <- a + 1
          b_used[a] <- b
          bootaux[[a]] <- thisfit
        } else {
          bootauxbad[[b]] <- thisfit
        }
        if (a == Brequired) break
      }

      print(paste0("Num of bootstrap sample sets used: ", b))
      if (a < Brequired) {
        warning(paste0("Only ", a, " bootstrap samples without variances
                   on constraint boundaries achieved, which was less
                   than the required number of ", Brequired))
        b_used <- b_used[1:a]
        Bstillneed <- Brequired - a
        badb_used <- setdiff(1:B, b_used)[1:Bstillneed]
        bootaux[(a + 1):Brequired] <- bootauxbad[badb_used]

      }
      varestmat <- t(sapply(1:Brequired, function(b) {
        bootaux[[b]]$var.est
      }, simplify = "array"))
    } else {
      b_used <- 1:B
      bootaux <- lapply(1:B, function(b) {
          thisfit <- alvm.fit(mainlm = bootobject[[b]],
              varselect = thevarselect,
              model = object$method, lambda = thelambda,
              nclust = object$hyperpar$nclustmethod,
              polypen = object$pentype, solver = object$solver,
              reduce2homosked = FALSE, ...)
      })
      varestmat <- t(sapply(1:B, function(b) {
        bootaux[[b]]$var.est
      }, simplify = "array"))
    }
  } else if (inherits(object, "anlvm.fit")) {
    if (rm_nonconverged) {

      bootaux <- vector("list", length = Brequired)
      a <- 0
      b_used <- rep(NA_integer_, Brequired)
      bootauxbad <- vector("list", length = B)

      for (b in 1:B) {
        thisfit <- anlvm.fit(mainlm = bootobject[[b]],
         varselect = thevarselect,
         g = object$fitinfo$g, cluster = (object$method == "cluster"),
         param.init = object$qlinfo$param.init,
         maxgridrows = object$qlinfo$maxgridrows,
         nconvstop = object$qlinfo$nconvstop,
         maxitql = object$qlinfo$maxitql,
         tolql = object$qlinfo$tolql, nestedql = object$qlinfo$nested,
         nclust = object$hyperpar$nclust, reduce2homosked = FALSE, ...)

        if (thisfit$qlinfo$converged) {
          a <- a + 1
          b_used[a] <- b
          bootaux[[a]] <- thisfit
        } else {
          bootauxbad[[b]] <- thisfit
        }
        if (a == Brequired) break
      }

      print(paste0("Num of bootstrap sample sets used: ", b))
      if (a < Brequired) {
        warning(paste0("Only ", a, " bootstrap samples resulted in
                   convergence of ANLVM estimation algorithm, which was less
                   than the required number of ", Brequired))
        b_used <- b_used[1:a]
        Bstillneed <- Brequired - a
        badb_used <- setdiff(1:B, b_used)[1:Bstillneed]
        bootaux[(a + 1):Brequired] <- bootauxbad[badb_used]
      }
      varestmat <- t(sapply(1:Brequired, function(b) {
        bootaux[[b]]$var.est
      }, simplify = "array"))
    } else {
      b_used <- 1:B
      bootaux <- lapply(1:B, function(b) {
        thisfit <- anlvm.fit(mainlm = bootobject[[b]],
            varselect = thevarselect,
            g = object$fitinfo$g, cluster = (object$method == "cluster"),
            param.init = object$qlinfo$param.init,
            maxgridrows = object$qlinfo$maxgridrows,
            nconvstop = object$qlinfo$nconvstop,
            maxitql = object$qlinfo$maxitql,
            tolql = object$qlinfo$tolql,
            nestedql = object$qlinfo$nested,
            nclust = object$hyperpar$nclust, reduce2homosked = FALSE, ...)
      })
      varestmat <- t(sapply(1:B, function(b) {
        bootaux[[b]]$var.est
      }, simplify = "array"))
    }
  }

  if (returnall) {
    value <- list("var.ests" = varestmat,
                  "avm" = bootaux, "b_used" = b_used)
  } else {
    value <- list("var.ests" = varestmat,
                  "avm" = NULL, "b_used" = b_used)
  }
  class(value) <- "bootavm"
  value
}

# Computes jackknife estimates of error variances,
#  to be used in BCa confidence intervals
jackavm <- function(object, retune = FALSE, ...) {

  if (!(inherits(object, "alvm.fit") ||
        inherits(object, "anlvm.fit")))
    stop("object should be of class alvm.fit or anlvm.fit")

  if (inherits(object, "alvm.fit")) {
    thelambda <- ifelse(retune, object$hyperpar$lambdamethod,
                        object$hyperpar$lambda)
  }

  X <- stats::model.matrix(object$ols)
  y <- stats::model.response(stats::model.frame(object$ols))
  n <- nrow(X)
  p <- ncol(X)

  olsjackknife <- lapply(1:n, function(j) {
    yleave1 <- y[-j]
    Xleave1 <- X[-j, , drop = FALSE]
    thelm <- stats::lm.fit(x = Xleave1, y = yleave1)
    list("y" = yleave1, "X" = Xleave1,
         "e" = stats::resid(thelm), "beta.hat" = stats::coef(thelm))
  })

  if (inherits(object, "alvm.fit")) {
    jackaux <- lapply(1:n, function(j) {
      alvm.fit(mainlm = olsjackknife[[j]],
               varselect = object$selectinfo$selectedcols,
               model = object$method, lambda = thelambda,
               nclust = object$hyperpar$nclust,
               polypen = object$pentype, solver = object$solver,
               reduce2homosked = FALSE, ...)
    })
  } else if (inherits(object, "anlvm.fit")) {
    jackaux <- lapply(1:n, function(j) {
      anlvm.fit(mainlm = olsjackknife[[j]],
              varselect = object$selectinfo$selectedcols,
              g = object$fitinfo$g, cluster = (object$method == "cluster"),
              param.init = object$qlinfo$param.init,
              maxgridrows = object$qlinfo$maxgridrows,
              nconvstop = object$qlinfo$nconvstop,
              maxitql = object$qlinfo$maxitql,
              tolql = object$qlinfo$tolql,
              nestedql = object$qlinfo$nested,
              nclust = object$hyperpar$nclust, reduce2homosked = FALSE, ...)
    })
  }

  if (inherits(object, "alvm.fit")) {
    fullvarests <- sapply(1:n, function(j) {
      drop(object$fitinfo$L %*% jackaux[[j]]$coef.est)
    }, simplify = "array")
  } else if (inherits(object, "anlvm.fit")) {
    fullvarests <- sapply(1:n, function(j) {
      vapply(1:n, function(k) {
        object$fitinfo$g(sum(object$fitinfo$Z[k, ] *
                               jackaux[[j]]$coef.est))
      }, 0)
    })
  }

  value <- list("jackvarests" = fullvarests, "jackavms" = jackaux)
  class(value) <- "jackavm"
  value
}

jackpoint <- function(jackobject, constol = 1e-10,
              aggfunc = c("mean", "median"), ...) {
  # Determine which jackknife fits have no estimates
  #  on constraint boundary
  n <- ncol(jackobject$jackvarests)
  aggfunc <- match.arg(aggfunc, c("mean", "median"))

  if (aggfunc == "mean") {
    aggfunc <- mean
  } else if (aggfunc == "median") {
    aggfunc <- stats::median
  }

  if (inherits(jackobject$jackavms[[1]], "alvm.fit")) {
    newvar.est <- vapply(1:n, function(i) {
      sdiff <- setdiff(which(jackobject$jackvarests[i, ] > 1.1 * constol), i)
      aggfunc(jackobject$jackvarests[i, sdiff])
    }, 0)
    newvar.est[is.nan(newvar.est)] <- constol

  } else if (inherits(jackobject$jackavms[[1]], "anlvm.fit")) {

    jackconverged <- which(vapply(1:n, function(j) {
      jackobject$jackavms[[j]]$qlinfo$converged
    }, TRUE))

    newvar.est <- vapply(1:n, function(i) {
      sdiff <- setdiff(jackconverged, i)
      aggfunc(jackobject$jackvarests[i, sdiff])
    }, 0)

  }
  unname(newvar.est)
}
