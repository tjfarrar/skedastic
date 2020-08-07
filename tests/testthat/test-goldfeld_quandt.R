context("goldfeld_quandt works for two lm examples across all argument permutations")


# theargs <- formals(goldfeld_quandt)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

test_that("parametric test: linear regression works with all combinations of formals", {
  theargs.par <- list("deflator" = c(NA, "speed", "crim", "2"), "method" = c("parametric"),
                  "prop_central" = c(1 / 3, 1 / 4), "group1prop" = c(1 / 2, 2 / 5),
                  "alternative" = c("greater", "less", "two.sided"),
                  "twosidedmethod" = c("doubled", "kulinskaya"),
                   "mainlm" = list(carslm, bostonlm))

  allargs.par <- expand.grid(theargs.par, stringsAsFactors = FALSE)
  allargs.par <- allargs.par[-which(vapply(1:nrow(allargs.par), function(i) allargs.par$deflator[i] == "speed" &
              !("speed" %in% colnames(model.matrix(allargs.par$mainlm[[i]]))), NA)), ]
  allargs.par <- allargs.par[-which(vapply(1:nrow(allargs.par), function(i) allargs.par$deflator[i] == "crim" &
              !("crim" %in% colnames(model.matrix(allargs.par$mainlm[[i]]))), NA)), ]

  pvals.par <- vapply(1:nrow(allargs.par), function(i) do.call(what = goldfeld_quandt,
              args = append(list("statonly" = FALSE),
              unlist(allargs.par[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals.par), function(i) expect_true(is.btwn01(pvals.par[i])))
})

ncars <- nrow(model.matrix(carslm))
nboston <- nrow(model.matrix(bostonlm))

test_that("nonparametric test with prob NA: linear regression works with all combinations of formals", {
  theargs.npar1 <- list("deflator" = c(NA, "speed", "crim", "2"), "method" = c("nonparametric"),
                      "restype" = c("ols", "blus"), "prob" = NA,
                      "alternative" = c("greater", "less", "two.sided"),
                      "twosidedmethod" = c("doubled", "kulinskaya"),
                      "mainlm" = list(carslm, bostonlm))

  allargs.npar1 <- expand.grid(theargs.npar1, stringsAsFactors = FALSE)
  allargs.npar1 <- allargs.npar1[-which(vapply(1:nrow(allargs.npar1),
              function(i) allargs.npar1$deflator[i] == "speed" &
            !("speed" %in% colnames(model.matrix(allargs.npar1$mainlm[[i]]))), NA)), ]
  allargs.npar1 <- allargs.npar1[-which(vapply(1:nrow(allargs.npar1),
              function(i) allargs.npar1$deflator[i] == "crim" &
            !("crim" %in% colnames(model.matrix(allargs.npar1$mainlm[[i]]))), NA)), ]

  pvals.npar1 <- vapply(1:nrow(allargs.npar1), function(i) do.call(what = goldfeld_quandt,
            args = append(list("statonly" = FALSE),
            unlist(allargs.npar1[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals.npar1), function(i) expect_true(is.btwn01(pvals.npar1[i])))
})

test_that("nonparametric test with prob not NA: linear regression works with all combinations of formals", {
  theargs.npar2 <- list("deflator" = c(NA, "speed", "crim", "2"), "method" = c("nonparametric"),
                        "restype" = c("ols", "blus"), "prob" = list(dpeakdat[[ncars]],
                    dpeakdat[[nboston]], dpeakdat[[ncars - 2]], dpeakdat[[nboston - 14]]),
                        "alternative" = c("greater", "less", "two.sided"),
                        "twosidedmethod" = c("doubled", "kulinskaya"),
                        "mainlm" = list(carslm, bostonlm))

  allargs.npar2 <- expand.grid(theargs.npar2, stringsAsFactors = FALSE)
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$deflator[i] == "speed" &
      !("speed" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$deflator[i] == "crim" &
      !("crim" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$restype[i] == "blus" &&
        length(allargs.npar2$prob[[i]]) != nboston - 14 &&
      !("speed" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$restype[i] == "ols" &&
        length(allargs.npar2$prob[[i]]) != nboston &&
      !("speed" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$restype[i] == "blus" &&
        length(allargs.npar2$prob[[i]]) != ncars - 2 &&
      !("crim" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]
  allargs.npar2 <- allargs.npar2[-which(vapply(1:nrow(allargs.npar2),
      function(i) allargs.npar2$restype[i] == "ols" &&
        length(allargs.npar2$prob[[i]]) != ncars &&
      !("crim" %in% colnames(model.matrix(allargs.npar2$mainlm[[i]]))), NA)), ]

  pvals.npar2 <- vapply(1:nrow(allargs.npar2), function(i) do.call(what = goldfeld_quandt,
            args = append(list("statonly" = FALSE), unlist(allargs.npar2[i, ],
                      recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals.npar2), function(i) expect_true(is.btwn01(pvals.npar2[i])))
})




