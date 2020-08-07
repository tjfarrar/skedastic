context("rackauskas_zuokas works for two lm examples across all argument permutations")

# theargs <- formals(rackauskas_zuokas)

test_that("linear regression works with all combinations of formals (data)", {
  carslm <- lm(dist ~ speed, data = cars)
  bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
  theargs <- list("alpha" = c(0, 15 / 32), "pvalmethod" = c("data"),
                  "R" = 2 ^ 14, "m" = 2 ^ 17, "sqZ" = c(TRUE, FALSE),
              "mainlm" = list(carslm, bostonlm))
  allargs <- expand.grid(theargs, stringsAsFactors = FALSE)

  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = rackauskas_zuokas,
              args = append(list("statonly" = FALSE, "seed" = 1234),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})

test_that("linear regression works with all combinations of formals (sim)", {
  skip_on_cran()
  carslm <- lm(dist ~ speed, data = cars)
  bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
  theargs <- list("alpha" = c(0, 15 / 32), "pvalmethod" = c("sim"),
                  "R" = 2 ^ 14, "m" = 2 ^ 4, "sqZ" = c(TRUE, FALSE),
                  "mainlm" = list(carslm, bostonlm))
  allargs2 <- expand.grid(theargs, stringsAsFactors = FALSE)

  pvals2 <- vapply(1:nrow(allargs2), function(i) do.call(what = rackauskas_zuokas,
            args = append(list("statonly" = FALSE, "seed" = 1234),
                          unlist(allargs2[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals2), function(i) expect_true(is.btwn01(pvals2[i])))
})
