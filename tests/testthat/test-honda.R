context("honda works for two lm examples across all argument permutations")

# theargs <- formals(honda)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

theargs1 <- list("deflator" = c(NA, "speed", "2"), "qfmethod" = c("imhof",
            "davies", "integrate"), "alternative" = c("two.sided", "greater", "less"),
            "twosidedmethod" = c("doubled", "kulinskaya"),
                "mainlm" = list(carslm))

theargs2 <- list("deflator" = c(NA, "crim", "2"), "qfmethod" = c("imhof",
            "davies", "integrate"), "alternative" = c("two.sided", "greater", "less"),
                 "twosidedmethod" = c("doubled", "kulinskaya"),
                 "mainlm" = list(bostonlm))


allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

test_that("carslm works with all combinations of formals", {
  pvals <- vapply(1:nrow(allargs1), function(i) do.call(what = honda,
              args = append(list("statonly" = FALSE),
              unlist(allargs1[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})

test_that("bostonlm works with all combinations of formals", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs2), function(i) do.call(what = honda,
        args = append(list("statonly" = FALSE),
        unlist(allargs2[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})
