context("horn works for two lm examples across all argument permutations")


# theargs <- formals(horn)

test_that("linear regression works with all combinations of formals", {
  skip_on_cran()
  carslm <- lm(dist ~ speed, data = cars)
  bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
  theargs <- list("deflator" = c(NA, "speed", "crim", "2"),
              "alternative" = c("two.sided",
              "greater", "less"), "restype" = c("ols", "blus"),
              "mainlm" = list(carslm, bostonlm))
  allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
  allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "speed" &
              !("speed" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
  allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "crim" &
              !("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]

  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = horn,
              args = append(list("statonly" = FALSE, "exact" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
