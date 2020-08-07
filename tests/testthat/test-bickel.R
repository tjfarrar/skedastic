context("bickel works for two lm examples across all argument permutations")


theargs <- formals(bickel)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
theargs <- list("fitmethod" = c("lm", "rlm"), "a" = c("identity", "2", "exp"),
                "b" = c("hubersq", "tanhsq"), "scale_invariant" = c(TRUE, FALSE),
                "k" = c(1.345, 1), "mainlm" = list(carslm, bostonlm))

allargs <- expand.grid(theargs, stringsAsFactors = FALSE)

test_that("linear regression works with all combinations of formals", {
  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = bickel,
                    args = append(list("mainlm" = carslm, "statonly" = FALSE),
                    unlist(allargs[i, ], recursive = TRUE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
