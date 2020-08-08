context("diblasi_bowman works for two lm examples across all argument permutations")


# ignorecov = TRUE (much faster!)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                 age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
theargs1 <- list("distmethod" = c("moment.match", "bootstrap"),
                "H" = list(0.08, seq(0.01, 0.13, 0.01)),
                ignorecov = c(TRUE),
                mainlm = list(carslm))
theargs2 <- list("distmethod" = c("moment.match", "bootstrap"),
                 "H" = list(0.08, seq(0.01, 0.13, 0.01)),
                 ignorecov = c(TRUE),
                 mainlm = list(bostonlm))

allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

allargs1 <- allargs1[-which(vapply(1:nrow(allargs1), function(i)
  (length(allargs1$H[[i]]) > 1), NA)), ]

test_that("diblasi-bowman works with all combinations of formals (carslm)", {
  pvals <- vapply(1:nrow(allargs1), function(i) do.call(what = diblasi_bowman,
      args = append(list("B" = 500L, "seed" = 1234,
      "statonly" = FALSE),
      unlist(allargs1[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})

test_that("diblasi-bowman works with all combinations of formals (bostonlm)", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs2), function(i) do.call(what = diblasi_bowman,
                args = append(list("B" = 500L, "seed" = 1234,
                                       "statonly" = FALSE),
                unlist(allargs2[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})

### ignorecov = FALSE

theargs3 <- list("distmethod" = c("moment.match", "bootstrap"),
                 "H" = list(0.08, seq(0.01, 0.13, 0.01)),
                 ignorecov = c(FALSE),
                 mainlm = list(carslm))
theargs4 <- list("distmethod" = c("moment.match", "bootstrap"),
                 "H" = list(0.08, seq(0.01, 0.13, 0.01)),
                 ignorecov = c(FALSE),
                 mainlm = list(bostonlm))

allargs3 <- expand.grid(theargs3, stringsAsFactors = FALSE)
allargs4 <- expand.grid(theargs4, stringsAsFactors = FALSE)

allargs3 <- allargs3[-which(vapply(1:nrow(allargs3), function(i)
  (length(allargs3$H[[i]]) > 1), NA)), ]

test_that("diblasi-bowman works with all combinations of formals (carslm)", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs3), function(i) do.call(what = diblasi_bowman,
            args = append(list("B" = 500L, "seed" = 1234, "statonly" = FALSE),
            unlist(allargs3[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})

test_that("diblasi-bowman works with all combinations of formals (bostonlm)", {
  skip_on_cran()
  pvals <- vapply(1:nrow(allargs4), function(i) do.call(what = diblasi_bowman,
            args = append(list("B" = 500L, "seed" = 1234, "statonly" = FALSE),
            unlist(allargs4[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
