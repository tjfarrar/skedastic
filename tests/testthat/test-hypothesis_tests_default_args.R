context("hypothesis tests work with default arguments for different kinds of mainlm argument")

htest <- c("anscombe", "bamset", "bickel",
           "breusch_pagan", "carapeto_holt", "cook_weisberg", "diblasi_bowman", "evans_king", "glejser",
           "goldfeld_quandt", "harrison_mccabe", "harvey", "honda", "horn",
           "li_yao", "rackauskas_zuokas", "simonoff_tsai",
           "szroeter", "verbyla", "white_lm", "wilcox_keselman", "zhou_etal")

test_that("simple linear regression works with default arguments and
          gives same p-value for lm or list mainlm argument", {
  skip_on_cran()
  carslm <- lm(dist ~ speed, data = cars)
  carslist <- list("y" = cars$dist, "X" = cbind(1, cars$speed))
  carslmpvals <- unlist(lapply(htest, function(h) do.call(what = h,
                  args = list("mainlm" = carslm))$p.value))
  carslistpvals <- unlist(lapply(htest, function(h) do.call(what = h,
                    args = list("mainlm" = carslist))$p.value))
  lapply(carslmpvals, function(x) expect_true(is.btwn01(x)))
  expect_equal(carslmpvals, carslistpvals)
})

test_that("simple linear regression through origin works with default
          arguments and gives same p-value for lm or list mainlm
          argument", {
            skip_on_cran()
            carslm0 <- lm(dist ~ 0 + speed, data = cars)
            carslist0 <- list("y" = cars$dist, "X" = as.matrix(cars$speed))
            carslm0pvals <- unlist(lapply(htest, function(h) do.call(what = h,
              args = list("mainlm" = carslm0))$p.value))
            carslist0pvals <- unlist(lapply(htest, function(h) do.call(what = h,
              args = list("mainlm" = carslist0))$p.value))
            lapply(carslm0pvals, function(x) expect_true(is.btwn01(x)))
            expect_equal(carslm0pvals, carslist0pvals)
})

test_that("lm (list) with NA's works with default arguments", {
  skip_on_cran()
  ChickWeight2 <- rbind(ChickWeight, NA)
  chicklistNA <- list("y" = ChickWeight2$weight, "X" = cbind(1,
                      ChickWeight2$Time))
  lapply(htest, function(h) expect_true(is.btwn01(do.call(what = h,
                      args = list("mainlm" = chicklistNA))$p.value)))
})

test_that("multiple linear regression works with default arguments and
          gives same p-value for lm or list mainlm argument", {
            skip_on_cran()
            chicklm <- lm(weight ~ Time + Diet, data = ChickWeight)
            chicklist <- list("y" = ChickWeight$weight,
                              "X" = stats::model.matrix(chicklm))
            chicklmpvals <- unlist(lapply(htest, function(h) do.call(what = h,
                      args = list("mainlm" = chicklm))$p.value))
            chicklistpvals <- unlist(lapply(htest, function(h) do.call(what = h,
                      args = list("mainlm" = chicklist))$p.value))
            lapply(chicklmpvals, function(x) expect_true(is.btwn01(x)))
            expect_equal(chicklmpvals, chicklistpvals)
})

htest_nosimonoff <- setdiff(htest, "simonoff_tsai")
test_that("large multiple linear regression (BostonHousing) works with default arguments", {
          skip_on_cran()
            bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                          age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
            bostonlmpvals <- unlist(lapply(htest_nosimonoff, function(h) do.call(what = h,
                                args = list("mainlm" = bostonlm))$p.value))
            lapply(bostonlmpvals, function(x) expect_true(is.btwn01(x)))
          })

test_that("multiple linear regression through origin works with default
          arguments and gives same p-value for lm or list mainlm
          argument", {
            skip_on_cran()
            chicklm0 <- lm(weight ~ 0 + Time + Diet, data = ChickWeight)
            chicklist0 <- list("y" = ChickWeight$weight,
                        "X" = stats::model.matrix(chicklm0))
            chicklm0pvals <- unlist(lapply(htest, function(h) do.call(what = h,
                      args = list("mainlm" = chicklm0))$p.value))
            chicklist0pvals <- unlist(lapply(htest, function(h) do.call(what = h,
                      args = list("mainlm" = chicklist0))$p.value))
            lapply(chicklm0pvals, function(x) expect_true(is.btwn01(x)))
            expect_equal(chicklm0pvals, chicklist0pvals)
})

test_that("multiple linear regression with random data in vectors rather than
           data frames works with default arguments and
           gives same p-value for lm or list mainlm argument", {
          skip_on_cran()
          myn <- 10
          myp <- 3
          myX <- matrix(data = runif(myn * myp),
                                 nrow = myn, ncol = myp)
          myy <- rnorm(myn, mean = rowMeans(myX), sd = rowMeans(myX))
          mylm <- lm(myy ~ myX)
          mylist <- list("y" = myy, "X" = cbind(1, myX))
          mylmpvals <- unlist(lapply(htest_nosimonoff, function(h) do.call(what = h,
                      args = list("mainlm" = mylm))$p.value))
          mylistpvals <- unlist(lapply(htest_nosimonoff, function(h) do.call(what = h,
                      args = list("mainlm" = mylist))$p.value))
          lapply(mylmpvals, function(x) expect_true(is.btwn01(x)))
          expect_equal(mylmpvals, mylistpvals)
})
