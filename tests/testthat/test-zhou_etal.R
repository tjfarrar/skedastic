context("zhou_etal works for two lm examples across all argument permutations")

theargs <- formals(zhou_etal)

carslm <- lm(dist ~ speed, data = cars)
ncars <- nrow(model.matrix(carslm))
carsauxmat <- cbind(1, matrix(data = runif(ncars * 5), nrow = ncars, ncol = 5))

bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
nboston <- nrow(model.matrix(bostonlm))
bostonauxmat <- cbind(1, matrix(data = runif(nboston * 5), nrow = nboston, ncol = 5))


theargs <- list("auxdesign" = list(NA, carsauxmat,
                bostonauxmat, "fitted.values"), "method" = c("pooled",
                "covariate-specific", "hybrid"), "Bperturbed" = 500L,
                "seed" = 1234,
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

test_that("linear regression works with all combinations of formals", {
  pvals <- unlist(lapply(1:nrow(allargs), function(i) do.call(what = zhou_etal,
                args = append(list("statonly" = FALSE),
                unlist(allargs[i, ], recursive = FALSE)))$p.value))
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
