context("simonoff_tsai works for two lm examples across all argument permutations")

theargs <- formals(simonoff_tsai)

carslm <- lm(dist ~ speed, data = cars)
ncars <- nrow(model.matrix(carslm))
carsauxmat <- cbind(1, matrix(data = runif(ncars * 5), nrow = ncars, ncol = 5))

bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
nboston <- nrow(model.matrix(bostonlm))
bostonauxmat <- cbind(1, matrix(data = runif(nboston * 5), nrow = nboston, ncol = 5))

theargs <- list("auxdesign" = list(NA, carsauxmat,
                bostonauxmat, "fitted.values"), "method" = c("mlr", "score"),
                "hetfun" = c("mult", "add", "logmult"),
                "basetest" = c("koenker", "cook_weisberg"),
                "bartlett" = c(TRUE, FALSE),
                "mainlm" = list(carslm, bostonlm))

allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
dim1 <- vapply(1:nrow(allargs), function(i) ifelse(is.null(dim(allargs$auxdesign[[i]])),
                0L, dim(allargs$auxdesign[[i]])[1]), NA_integer_)
allargs <- allargs[-which(vapply(1:nrow(allargs),
              function(i) dim1[i] == nboston &
              "speed" %in% colnames(model.matrix(allargs$mainlm[[i]])), NA)), ]
dim1 <- vapply(1:nrow(allargs), function(i) ifelse(is.null(dim(allargs$auxdesign[[i]])),
                0L, dim(allargs$auxdesign[[i]])[1]), NA_integer_)
allargs <- allargs[-which(vapply(1:nrow(allargs),
              function(i) dim1[i] == ncars &
              "crim" %in% colnames(model.matrix(allargs$mainlm[[i]])), NA)), ]
nmodel <- vapply(allargs$mainlm, function(ell) nrow(model.matrix(ell)), NA_integer_)

# which1 <- which(allargs$method == "mlr" & nmodel == 506)
which2 <- which(allargs$hetfun == "logmult" & allargs$auxdesign == "fitted.values")
which3 <- which(allargs$method == "score" & allargs$hetfun == "add")
which4 <- which(allargs$method == "score" & allargs$bartlett)
which5 <- which(allargs$hetfun == "add" & allargs$bartlett)
which6 <- which(allargs$hetfun == "logmult" & nmodel == 506 & is.na(allargs$auxdesign))
which7 <- which(allargs$method == "mlr" & allargs$basetest == "koenker")
allargs <- allargs[-unique(c(which2, which3, which4, which5, which6,
                             which7)), ]


test_that("linear regression works with all combinations of formals", {
  pvals <- vapply(1:23, function(i)
    suppressWarnings(do.call(what = simonoff_tsai,
                args = append(list("statonly" = FALSE),
                unlist(allargs[i, ], recursive = FALSE)))$p.value), NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})

test_that("linear regression works with all combinations of formals", {
  skip_on_cran()
  pvals <- vapply(24:nrow(allargs), function(i)
    suppressWarnings(do.call(what = simonoff_tsai,
                args = append(list("statonly" = FALSE),
                unlist(allargs[i, ], recursive = FALSE)))$p.value), NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
