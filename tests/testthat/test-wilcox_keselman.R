context("wilcox_keselman works for two lm examples across all argument permutations")

# theargs <- formals(wilcox_keselman)

test_that("linear regression works with all combinations of formals", {
  carslm <- lm(dist ~ speed, data = cars)
  bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
  theargs <- list("gammapar" = c(0.2, 0.4), "p.adjust.method" = c("none", "holm"),
                  "seed" = c(NA, 2), "matchWRS" = c(TRUE, FALSE),
              "mainlm" = list(carslm, bostonlm))
  allargs <- expand.grid(theargs, stringsAsFactors = FALSE)

  pvals <- unlist(lapply(1:nrow(allargs), function(i) do.call(what = wilcox_keselman,
              args = append(list("statonly" = FALSE, "B" = 500L, "rqwarn" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value))
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
