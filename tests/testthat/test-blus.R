context("blus works for cars lm and boston lm across all argument permutations")

# theargs <- formals(blus)

carslm <- lm(dist ~ speed, data = cars)
n <- nrow(model.matrix(carslm))
p <- ncol(model.matrix(carslm))
theargs1 <- list("omit" = c("first", "last", "random"), "exhaust" = c(NA, 3),
                 "keepNA" = c(TRUE, FALSE), "seed" = c(1234, NA))
theargs2 <- list("omit" = list(c(2, 20)), "exhaust" = c(NA, 3),
                 "keepNA" = c(TRUE, FALSE), "seed" = c(1234, NA))

allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

test_that("simple linear regression (carslm) works with all combinations of formals", {
  carsblus1 <- lapply(1:nrow(allargs1), function(i) do.call(what = blus,
                    args = append(list("mainlm" = carslm), allargs1[i, ])))
  lapply(1:length(carsblus1), function(i)
    expect_true(sum(is.finite(carsblus1[[i]])) == n - p))

  carsblus2 <- lapply(1:nrow(allargs2), function(i) do.call(what = blus,
           args = append(list("mainlm" = carslm), unlist(allargs2[i, ],
                                                      recursive = FALSE))))
  lapply(1:length(carsblus2), function(i)
    expect_true(sum(is.finite(carsblus2[[i]])) == n - p))
})


bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
n <- nrow(model.matrix(bostonlm))
p <- ncol(model.matrix(bostonlm))
theargs1 <- list("omit" = c("first", "last", "random"), "exhaust" = c(NA, 5),
                 "keepNA" = c(TRUE, FALSE), "seed" = c(1234, NA))
theargs2 <- list("omit" = list(seq(from = 5, by = 2, length.out = 14)), "exhaust" = c(NA, 5),
                 "keepNA" = c(TRUE, FALSE), "seed" = c(1234, NA))

allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

test_that("simple linear regression (bostonlm) works with all combinations of formals", {
  bostonblus1 <- lapply(1:nrow(allargs1), function(i) do.call(what = blus,
        args = append(list("mainlm" = bostonlm), allargs1[i, ])))
  lapply(1:length(bostonblus1), function(i)
    expect_true(sum(is.finite(bostonblus1[[i]])) == n - p))

  bostonblus2 <- lapply(1:nrow(allargs2), function(i) do.call(what = blus,
        args = append(list("mainlm" = bostonlm), unlist(allargs2[i, ], recursive = FALSE))))
  lapply(1:length(bostonblus2), function(i)
    expect_true(sum(is.finite(bostonblus2[[i]])) == n - p))
})
