context("hypothesis tests work with default arguments for different kinds of mainlm argument")

htest <- c("anscombe", "bamset", "bickel",
           "breusch_pagan", "carapeto_holt", "cook_weisberg", "diblasi_bowman", "evans_king", "glejser",
           "goldfeld_quandt", "harvey", "honda", "horn",
           "li_yao", "rackauskas_zuokas", "simonoff_tsai",
           "szroeter", "verbyla", "white_lm", "wilcox_keselman", "zhou_etal")

test_that("simple lm works with default arguments", {

  carslm <- lm(dist ~ speed, data = cars)
  lapply(htest, function(h) expect_true(is.finite(do.call(what = h,
                          args = list("mainlm" = carslm))$p.value)))
})

test_that("simple lm thru origin works with default arguments", {
  carslm0 <- lm(dist ~ 0 + speed, data = cars)
  lapply(htest, function(h) expect_true(is.finite(do.call(what = h,
        args = list("mainlm" = carslm0))$p.value)))
})

test_that("simple lm (list) works with default arguments", {
  carslist <- list("y" = cars$dist, "X" = cbind(1, cars$speed))
  lapply(htest, function(h) expect_true(is.finite(do.call(what = h,
              args = list("mainlm" = carslist))$p.value)))
})

test_that("simple lm (list) thru origin works with default arguments", {
  carslist0 <- list("y" = cars$dist, "X" = as.matrix(cars$speed))
  lapply(htest, function(h) expect_true(is.finite(do.call(what = h,
              args = list("mainlm" = carslist0))$p.value)))
})

# May need to move to separate file and skip_on_cran() because time-consuming
# test_that("lm (list) with NA's works with default arguments", {
#   ChickWeight2 <- rbind(ChickWeight, NA)
#   chicklistNA <- list("y" = ChickWeight2$weight, "X" = cbind(1,
#                                     ChickWeight2$Time))
#   lapply(htest, function(h) expect_true(is.finite(do.call(what = h,
#         args = list("mainlm" = chicklistNA))$p.value)))
# })

# ChickWeight2 <- rbind(ChickWeight, NA)
# themainlm <- vector("list", 8)
# theresults <- vector("list", 8)
# themainlm[[1]] <- lm(weight ~ Time, data = ChickWeight)
# themainlm[[2]] <- list("y" = ChickWeight$weight, "X" = cbind(1, ChickWeight$Time))
# themainlm[[3]] <- list("y" = ChickWeight$weight, "X" = cbind(1, ChickWeight$Time),
#                        "e" = themainlm[[1]]$residuals)
# themainlm[[4]] <- lm(weight ~ Time + Diet, data = ChickWeight)
# themainlm[[5]] <- lm(weight ~ 0 + Time, data = ChickWeight)
# myn <- 10
# myp <- 4
# myX <- cbind(1, matrix(data = runif(myn * myp), nrow = myn, ncol = myp))
# myy <- rnorm(myn, mean = rowMeans(myX), sd = rowMeans(myX))
# themainlm[[6]] <- lm(myy ~ myX)
# themainlm[[7]] <- list("y" = myy, "X" = myX)
# themainlm[[8]] <- lm(weight ~ Time, data = ChickWeight2)
#
