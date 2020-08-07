context("Dtrend functions work in a variety of settings")

# dDtrend

test_that("dDtrend: test some unusual sequences", {
  expect_error(dDtrend(k = NULL, n = 8))
  expect_error(dDtrend(k = "all", n = 12))
  expect_error(dDtrend(k = "alll", n = 5))
  expect_error(dDtrend(k = 1.5, n = 5))
  expect_error(dDtrend(k = 1, n = 5.5))
  expect_equal(dDtrend(k = "all", n = 4)[c(2, 4)],
               dDtrend(k = c(2, 6), n = 4))
})

# pDtrend

test_that("dDtrend: test some unusual sequences", {
  expect_error(pDtrend(k = NULL, n = 8))
  expect_error(pDtrend(k = "all", n = 12))
  expect_error(pDtrend(k = "all", n = 5, exact = FALSE))
  expect_error(pDtrend(k = 1, n = 5.5))
  expect_equal(pDtrend(k = "all", n = 4)[c(2, 4)],
               pDtrend(k = c(2, 6), n = 4))
})
