context("anlvm.fit works for cars lm and boston lm across several argument permutations")

set.seed(234)
carslm <- lm(dist ~ speed, data = cars)


theargs1 <- list("g" = c("sq", "exp"),
                 "varselect" = c("none", "hettest",
                      "cv.linear", "cv.cluster",
                      "qgcv.linear", "qgcv.cluster"))

theargs2 <- list("g" = "sq",
                 "cluster" = TRUE,
                 "nclust" = c("elbow.swd", "elbow.mwd", "elbow.both"),
                 "varselect" = c("none", "hettest",
                                 "cv.linear", "cv.cluster",
                                 "qgcv.linear", "qgcv.cluster"))

allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

test_that("simple linear regression (carslm) works with all combinations of formals", {
  carsalvm1 <- lapply(1:nrow(allargs1), function(i) {
    # print(i)
    do.call(what = anlvm.fit,
                    args = append(list("mainlm" = carslm), allargs1[i, ]))
  })
  lapply(1:length(carsalvm1), function(i)
    expect_true(all(carsalvm1[[i]]$var.est > 0)))

  carsalvm2 <- lapply(1:nrow(allargs2), function(i) {
    # print(i)
    do.call(what = anlvm.fit,
                    args = append(list("mainlm" = carslm), allargs2[i, ]))
  })
  lapply(1:length(carsalvm2), function(i)
    expect_true(all(carsalvm2[[i]]$var.est > 0)))

})


bostonlm <- lm(cmedv ~ crim + zn + indus + chas + nox + rm +
    age + dis + rad + tax + ptratio + b + lstat,
    data = BostonHousing2)

theargs1 <- list("g" = c("sq", "exp"),
                 "varselect" = c("none", "hettest",
                                 "qgcv.linear"),
                 "hettest" = "breusch_pagan")

theargs2 <- list("g" = "sq",
                 "cluster" = TRUE,
                 "nclust" = c("elbow.swd", "elbow.mwd", "elbow.both"),
                 "varselect" = c("none", "hettest",
                                 "qgcv.linear"),
                 "hettest" = "breusch_pagan")

allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)


test_that("simple linear regression (bostonlm) works with all combinations of formals", {
  skip_on_cran()
  bostonalvm1 <- lapply(1:nrow(allargs1), function(i) {
    # print(i)
    do.call(what = anlvm.fit,
                args = append(list("mainlm" = bostonlm), allargs1[i, ]))
  })
  lapply(1:length(bostonalvm1), function(i)
    expect_true(all(bostonalvm1[[i]]$var.est > 0)))

  bostonalvm2 <- lapply(1:nrow(allargs2), function(i) {
    # print(i)
    do.call(what = anlvm.fit,
      args = append(list("mainlm" = bostonlm), allargs2[i, ]))
    })
  lapply(1:length(bostonalvm2), function(i)
    expect_true(all(bostonalvm2[[i]]$var.est > 0)))

})
