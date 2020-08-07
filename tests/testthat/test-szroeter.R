context("szroeter works for two lm examples across all argument permutations")

# theargs <- formals(szroeter)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
theargs <- list("deflator" = c(NA, "speed", "crim", "2"),
                "h" = list(SKH), "qfmethod" = c("imhof", "davies", "integrate"),
                "mainlm" = list(carslm, bostonlm))
allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "speed" &
            !("speed" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "crim" &
            !("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]


test_that("linear regression works with all combinations of formals", {
  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = szroeter,
              args = append(list("statonly" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
