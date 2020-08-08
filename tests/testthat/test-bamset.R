context("bamset works for two lm examples across all argument permutations")

theargs <- formals(bamset)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

theargs <- list("k" = c(2, 5), "deflator" = c(NA, "speed", "crim", "2"),
                "correct" = c(TRUE, FALSE), "omitatmargins" = c(TRUE, FALSE),
                "omit" = c("first", "last", "random", NA),
                mainlm = list(carslm, bostonlm))

allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
allargs <- allargs[-which(is.na(allargs$omit) & allargs$omitatmargins == FALSE), ]
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "speed" &
              !("speed" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "crim" &
              !("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]

test_that("linear regression works with all combinations of formals", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = bamset,
                  args = append(list("statonly" = FALSE),
                  unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
