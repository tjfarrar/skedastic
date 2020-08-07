context("hypothesis tests that require only OLS residuals work when mainlm is a
        list object containing only the OLS residuals")

etest <- c("li_yao", "li_yao", "rackauskas_zuokas", "rackauskas_zuokas")

test_that("simple linear regression works when mainlm argument consists
          only of OLS residuals", {
  carslm <- lm(dist ~ speed, data = cars)
  carslist <- list("e" = carslm$residuals)
  eargs <- list(list("mainlm" = carslist, "method" = "alrt"),
                list("mainlm" = carslist, "method" = "cvt", "baipanyin" = FALSE),
                list("mainlm" = carslist, "alpha" = 0),
                list("mainlm" = carslist, "alpha" = 15 / 32))
  carslistpvals <- vapply(1:length(etest), function(i) do.call(what = etest[i],
                    args = eargs[[i]])$p.value, NA_real_)
  lapply(carslistpvals, function(x) expect_true(is.btwn01(x)))
})
