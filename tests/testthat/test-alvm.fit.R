context("alvm.fit works for cars lm and boston lm across several argument permutations")

carslm <- lm(dist ~ speed, data = cars)
theargs1 <- list("model" = "linear",
                 "varselect" = c("none", "hettest",
                                 "cv.linear", "cv.cluster",
                                 "qgcv.linear", "qgcv.cluster"))
theargs2 <- list("model" = "cluster",
                 "nclust" = c("elbow.swd", "elbow.mwd", "elbow.both",
                              "foldcv"),
                 "varselect" = c("none", "hettest",
                                 "cv.linear", "cv.cluster",
                                 "qgcv.linear", "qgcv.cluster"))
theargs3 <- list("model" = "polynomial",
                 "lambda" = c("foldcv", "qgcv"),
                 "polypen" = c("L2", "L1"),
                 "d" = c(1, 2, 3),
                 "cvoption" = c("testsetols", "partitionres"))
theargs4 <- list("model" = "spline",
                 "lambda" = c("foldcv", "qgcv"),
                 "cvoption" = c("testsetols", "partitionres"))
theargs5 <- list("model" = c("basic", "homoskedastic"))


allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)
allargs3 <- expand.grid(theargs3, stringsAsFactors = FALSE)
allargs4 <- expand.grid(theargs4, stringsAsFactors = FALSE)
allargs5 <- expand.grid(theargs5, stringsAsFactors = FALSE)

test_that("simple linear regression (carslm) works with all combinations of formals", {
  skip_on_cran()
  carsalvm1 <- lapply(1:nrow(allargs1), function(i) do.call(what = alvm.fit,
                    args = append(list("mainlm" = carslm), allargs1[i, ])))
  lapply(1:length(carsalvm1), function(i)
    expect_true(all(carsalvm1[[i]]$var.est > 0)))

  carsalvm2 <- lapply(1:nrow(allargs2), function(i) do.call(what = alvm.fit,
                    args = append(list("mainlm" = carslm), allargs2[i, ])))
  lapply(1:length(carsalvm2), function(i)
    expect_true(all(carsalvm2[[i]]$var.est > 0)))

  carsalvm3 <- lapply(1:nrow(allargs3), function(i) {
    # print(i)
    do.call(what = alvm.fit,
        args = append(list("mainlm" = carslm), allargs3[i, ]))
    })
  lapply(1:length(carsalvm3), function(i)
    expect_true(all(carsalvm3[[i]]$var.est > 0)))

  carsalvm4 <- lapply(1:nrow(allargs4), function(i) {
    # print(i)
    do.call(what = alvm.fit,
        args = append(list("mainlm" = carslm), allargs4[i, ]))
    })
  lapply(1:length(carsalvm4), function(i)
    expect_true(all(carsalvm4[[i]]$var.est > 0)))

  carsalvm5 <- lapply(1:nrow(allargs5), function(i) {
    do.call(what = alvm.fit,
        args = append(list("mainlm" = carslm), allargs5[i, , drop = FALSE]))
    })
  lapply(1:length(carsalvm5), function(i)
    expect_true(all(carsalvm5[[i]]$var.est > 0)))

})


bostonlm <- lm(cmedv ~ crim + zn + indus + chas + nox + rm +
    age + dis + rad + tax + ptratio + b + lstat,
    data = BostonHousing2)

theargs1 <- list("model" = "linear",
                 "varselect" = c("none", "hettest",
                                 "qgcv.linear"))
theargs2 <- list("model" = "cluster",
                 "nclust" = c("elbow.swd", "elbow.mwd", "elbow.both",
                              "foldcv"),
                 "varselect" = c("none", "hettest",
                                 "qgcv.linear"))
theargs3 <- list("model" = "polynomial",
                 "lambda" = c("foldcv", "qgcv"),
                 "polypen" = c("L2", "L1"),
                 "d" = c(1L, 2L),
                 "cvoption" = c("testsetols", "partitionres"))
theargs4 <- list("model" = c("basic", "homoskedastic"))


allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)
allargs3 <- expand.grid(theargs3, stringsAsFactors = FALSE)
allargs4 <- expand.grid(theargs4, stringsAsFactors = FALSE)

test_that("simple linear regression (bostonlm) works with all combinations of formals", {
  skip_on_cran()
  bostonalvm1 <- lapply(1:nrow(allargs1), function(i) do.call(what = alvm.fit,
                args = append(list("mainlm" = bostonlm), allargs1[i, ])))
  lapply(1:length(bostonalvm1), function(i)
    expect_true(all(bostonalvm1[[i]]$var.est > 0)))

  bostonalvm2 <- lapply(1:nrow(allargs2), function(i) do.call(what = alvm.fit,
                args = append(list("mainlm" = bostonlm), allargs2[i, ])))
  lapply(1:length(bostonalvm2), function(i)
    expect_true(all(bostonalvm2[[i]]$var.est > 0)))

  bostonalvm3 <- lapply(1:nrow(allargs3), function(i) {
    # print(i)
    do.call(what = alvm.fit,
          args = append(list("mainlm" = bostonlm), allargs3[i, ]))
    })
  lapply(1:length(bostonalvm3), function(i)
    expect_true(all(bostonalvm3[[i]]$var.est > 0)))

  bostonalvm4 <- lapply(1:nrow(allargs4), function(i) do.call(what = alvm.fit,
                args = append(list("mainlm" = bostonlm),
                              allargs4[i, , drop = FALSE])))
  lapply(1:length(bostonalvm4), function(i)
    expect_true(all(bostonalvm4[[i]]$var.est > 0)))
})
