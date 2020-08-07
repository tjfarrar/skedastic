context("peaks functions work in a variety of settings")

# countpeaks

carslm <- lm(dist ~ speed, data = cars)

test_that("countpeaks: test some unusual sequences", {
  expect_error(countpeaks(rep(NA_real_, 100)))
  expect_true(countpeaks(c(rep(NaN, 50), 1)) == 0L)
  expect_true(countpeaks(1) == 0L)
  expect_error(countpeaks(double(0)))
  expect_true(countpeaks(c(-Inf, 3, Inf)) == 2L)
  expect_true(countpeaks(blus(mainlm = carslm)) > 0L)
})

# dpeak

test_that("dpeak: test correctness of data, n <= 170", {
  dpeaknodat <- dpeak(0:169, 170, usedata = FALSE)
  expect_true(all.equal(dpeaknodat, dpeakdat[[170]]))

  ind <- sample(169, 30, replace = TRUE)
  dpeaknodat2 <- dpeak(ind, 170, usedata = FALSE)
  expect_true(all.equal(dpeaknodat2, dpeakdat[[170]][ind + 1]))

})

test_that("dpeak: test correctness of data, n > 170", {
  skip_on_cran()
  dpeaknodat <- dpeak(c(2, 6), 180, usedata = FALSE)
  expect_equal(dpeaknodat, dpeakdat[[180]][c(2, 6) + 1])
})

test_that("dpeak: errors expected", {
  expect_error(dpeak(k = 3.5, n = 4))
  expect_error(dpeak(k = 0:5, n = 5))
  expect_error(dpeak(k = 0, n = 1001, usedata = TRUE))
})

# ppeak

test_that("ppeak: correctness, n <= 170, lower tail", {
  ppeaknodat <- ppeak(k = c(1, 3, 5, 7), n = 150, usedata = FALSE,
                      lower.tail = TRUE)
  ppeakdat <- ppeak(k = c(1, 3, 5, 7), n = 150, usedata = TRUE,
                    lower.tail = TRUE)
  noppeakdat <- vapply(c(1, 3, 5, 7), function(j)
    sum(dpeakdat[[150]][1:(j + 1)]), NA_real_)
  expect_true(all.equal(ppeaknodat, noppeakdat))
  expect_true(all.equal(ppeakdat, noppeakdat))
})

test_that("ppeak: correctness, n <= 170, upper tail", {
  ppeaknodat <- ppeak(k = c(8, 5, 3), n = 67, usedata = FALSE,
                      lower.tail = FALSE)
  ppeakdat <- ppeak(k = c(8, 5, 3), n = 67, usedata = TRUE,
                    lower.tail = FALSE)
  noppeakdat <- vapply(c(8, 5, 3), function(j)
    sum(dpeakdat[[67]][(j + 1):66]), NA_real_)
  expect_true(all.equal(ppeaknodat, noppeakdat))
  expect_true(all.equal(ppeakdat, noppeakdat))
})

test_that("ppeak: correctness, n > 170, lower tail", {
  skip_on_cran()
  ppeaknodat <- ppeak(k = c(2, 6), n = 175, usedata = FALSE,
                      lower.tail = TRUE)
  ppeakdat <- ppeak(k = c(2, 6), n = 175, usedata = TRUE,
                    lower.tail = TRUE)
  noppeakdat <- vapply(c(2, 6), function(j)
    sum(dpeakdat[[175]][1:(j + 1)]), NA_real_)
  expect_true(all.equal(ppeaknodat, noppeakdat))
  expect_true(all.equal(ppeakdat, noppeakdat))
})

test_that("ppeak: correctness, n > 170, upper tail", {
  skip_on_cran()
  ppeaknodat <- ppeak(k = c(19, 4), n = 283, usedata = FALSE,
                      lower.tail = FALSE)
  ppeakdat <- ppeak(k = c(19, 4), n = 283, usedata = TRUE,
                    lower.tail = FALSE)
  noppeakdat <- vapply(c(19, 4), function(j)
    sum(dpeakdat[[283]][(j + 1):282]), NA_real_)
  expect_true(all.equal(ppeaknodat, noppeakdat))
  expect_true(all.equal(ppeakdat, noppeakdat))
})
