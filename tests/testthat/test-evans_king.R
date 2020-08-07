context("evans_king works for two lm examples across all argument permutations")

# theargs <- formals(evans_king)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

theargs_GLS1 <- list("deflator" = c(NA, "speed", "2"), "qfmethod" = c("imhof",
                "davies", "integrate"), "method" = c("GLS"),
                "lambda_star" = c(2.5, 5), "mainlm" = list(carslm))

theargs_GLS2 <- list("deflator" = c(NA, "crim", "2"), "qfmethod" = c("imhof",
                "davies", "integrate"), "method" = c("GLS"),
                "lambda_star" = c(2.5, 5), "mainlm" = list(bostonlm))

theargs_LM1 <- list("deflator" = c(NA, "speed", "2"), "qfmethod" = c("imhof",
                "davies", "integrate"), "method" = c("LM"),
                "lambda_star" = c(2.5, 5), "mainlm" = list(carslm))

theargs_LM2 <- list("deflator" = c(NA, "crim", "2"), "qfmethod" = c("imhof",
                "davies", "integrate"), "method" = c("LM"),
                "lambda_star" = c(2.5, 5), "mainlm" = list(bostonlm))

allargs_GLS1 <- expand.grid(theargs_GLS1, stringsAsFactors = FALSE)
allargs_GLS2 <- expand.grid(theargs_GLS2, stringsAsFactors = FALSE)

allargs_LM1 <- expand.grid(theargs_LM1, stringsAsFactors = FALSE)
allargs_LM2 <- expand.grid(theargs_LM2, stringsAsFactors = FALSE)

test_that("evans_king GLS method works with all combinations of formals; carslm", {
  pvals <- vapply(1:nrow(allargs_GLS1), function(i) do.call(what = evans_king,
      args = append(list("statonly" = FALSE),
      unlist(allargs_GLS1[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})


test_that("evans_king GLS method works with all combinations of formals; bostonlm", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs_GLS2), function(i) do.call(what = evans_king,
              args = append(list("statonly" = FALSE),
              unlist(allargs_GLS2[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})

test_that("evans_king LM1 method works with all combinations of formals; carslm", {
  pvals <- vapply(1:nrow(allargs_LM1), function(i) do.call(what = evans_king,
  args = append(list("statonly" = FALSE),
  unlist(allargs_LM1[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})

test_that("evans_king LM2 method works with all combinations of formals; carslm", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs_LM1), function(i) do.call(what = evans_king,
  args = append(list("statonly" = FALSE),
  unlist(allargs_LM2[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})
