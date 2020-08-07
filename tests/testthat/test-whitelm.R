context("white_lm works for two lm examples across all argument permutations")

theargs <- formals(white_lm)

carslm <- lm(dist ~ speed, data = cars)
ncars <- nrow(model.matrix(carslm))
carsauxmat <- cbind(1, matrix(data = runif(ncars * 5), nrow = ncars, ncol = 5))

bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
nboston <- nrow(model.matrix(bostonlm))
bostonauxmat <- cbind(1, matrix(data = runif(nboston * 5), nrow = nboston, ncol = 5))

theargs <- list("interactions" = c(TRUE, FALSE),
                "mainlm" = list(carslm, bostonlm))

allargs <- expand.grid(theargs, stringsAsFactors = FALSE)

test_that("linear regression works with all combinations of formals", {
  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = white_lm,
                args = append(list("statonly" = FALSE),
                unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
