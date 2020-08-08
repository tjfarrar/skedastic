context("carapeto_holt works for two lm examples across all argument permutations")


# theargs <- formals(carapeto_holt)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
theargs <- list("deflator" = c(NA, "speed", "crim", "2"), "prop_central" = c(1 / 2, 1 / 3,
                  1 / 4), "group1prop" = c(1 / 2, 1 / 3),
                "qfmethod" = c("imhof", "davies", "integrate"),
                "alternative" = c("greater", "less", "two.sided"),
                "twosidedmethod" = c("doubled", "kulinskaya"),
                "mainlm" = list(carslm, bostonlm))

allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "speed" &
                                   !("speed" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$deflator[i] == "crim" &
                                   !("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
# allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$qfmethod[i] == "integrate" &
#                                    ("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))), NA)), ]
# allargs <- allargs[-which(vapply(1:nrow(allargs), function(i) allargs$qfmethod[i] == "davies" &
#                                    ("crim" %in% colnames(model.matrix(allargs$mainlm[[i]]))) &
#                                    allargs$alternative[i] == "two.sided" &
#                                    allargs$twosidedmethod[i] == "kulinskaya", NA)), ]


test_that("linear regression works with all combinations of formals", {
  skip_on_cran()
  pvals <- vapply(1:324, function(i) do.call(what = carapeto_holt,
              args = append(list("statonly" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  unlist(lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i]))))
})

test_that("linear regression works with all combinations of formals", {
  skip_on_cran()
  pvals <- vapply(325:nrow(allargs), function(i) do.call(what = carapeto_holt,
              args = append(list("statonly" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  unlist(lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i]))))
})
