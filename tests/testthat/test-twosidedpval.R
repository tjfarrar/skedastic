context("twosidedpval works for an assortment of discrete and continuous distributions")

### t.test (symmetrical, so methods should give same result)

n <- c(1200, 1000)
myttest <- vector("list", 2)
normdat <- list("r0" = list(rnorm(n = n[1]), rnorm(n = n[2])),
                "r1" = list(rnorm(n = n[1], mean = 0.2), rnorm(n = n[2])),
                "r2" = list(rnorm(n = n[1], sd = 1.25), rnorm(n = n[2])))
myttest[[1]] <- stats::t.test(x = normdat$r0[[1]], y = normdat$r0[[2]],
                         alternative = "two.sided", var.equal = TRUE)
myttest[[2]] <- stats::t.test(x = normdat$r1[[1]], y = normdat$r1[[2]],
                         alternative = "two.sided", var.equal = TRUE)

test_that("t-tests give same value for `doubled` and `kulinskaya` methods", {
  skip_on_cran()
  pvals <- vapply(c("doubled", "kulinskaya"), function(m)
    twosidedpval(q = abs(myttest[[1]]$statistic), CDF = pt, continuous = TRUE,
    method = m, Aloc = 0, df = n[1] + n[2] - 2), NA_real_, USE.NAMES = FALSE)
  lapply(pvals, function(p) expect_equal(p, myttest[[1]]$p.value))
  lapply(pvals, function(p) expect_true(is.btwn01(p)))

  pvals <- vapply(c("doubled", "kulinskaya"), function(m)
    twosidedpval(q = abs(myttest[[2]]$statistic), CDF = pt, continuous = TRUE,
    method = m, Aloc = 0, df = n[1] + n[2] - 2), NA_real_, USE.NAMES = FALSE)
  lapply(pvals, function(p) expect_equal(p, myttest[[2]]$p.value))
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})


### F test of variance
myFtest <- vector("list", 2)
myFtest[[1]] <- stats::var.test(x = normdat$r0[[1]], y = normdat$r0[[2]],
                            alternative = "two.sided")
myFtest[[2]] <- stats::var.test(x = normdat$r2[[1]], y = normdat$r2[[2]],
                                alternative = "two.sided")

test_that("F-tests give different values for `doubled` and `kulinskaya` methods", {
  pvals <- vapply(c("doubled", "kulinskaya"), function(m)
    twosidedpval(q = myFtest[[1]]$statistic, CDF = pf, continuous = TRUE,
                 method = m, Aloc = 1, df1 = n[1] - 1, df2 = n[2] - 1),
                 NA_real_, USE.NAMES = FALSE)
  expect_equal(pvals[1], myFtest[[1]]$p.value)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))

  pvals <- vapply(c("doubled", "kulinskaya"), function(m)
    twosidedpval(q = myFtest[[2]]$statistic, CDF = pf, continuous = TRUE,
                 method = m, Aloc = 0, df1 = n[1] - 1, df2 = n[2] - 1),
                 NA_real_, USE.NAMES = FALSE)
  expect_equal(pvals[1], myFtest[[2]]$p.value)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})

### Binomial test

binomdat <- list("r0" = rbinom(n = 1, size = n[1], prob = 0.4),
                 "r1" = rbinom(n = 1, size = n[1], prob = 0.3))
mybtest <- vector("list", 2)
mybtest[[1]] <- stats::binom.test(x = binomdat$r0, n = n[1], p = 0.4)
mybtest[[2]] <- stats::binom.test(x = binomdat$r1, n = n[1], p = 0.4)

test_that("Binomial tests give different values for `doubled`, `kulinskaya`, and `minlikelihood` methods", {
  pvals <- vapply(c("doubled", "kulinskaya", "minlikelihood"), function(m)
    twosidedpval(q = mybtest[[1]]$statistic, CDF = pbinom, continuous = FALSE,
                 method = m, Aloc = n[1] * 0.4, supportlim = c(0, n[1]),
                 size = n[1], prob = 0.4), NA_real_, USE.NAMES = FALSE)
  expect_equal(pvals[3], mybtest[[1]]$p.value)
  pvalsmin_nosupport <- twosidedpval(q = mybtest[[1]]$statistic, CDF = pbinom,
                        continuous = FALSE, method = "m", Aloc = n[1] * 0.4,
                        size = n[1], prob = 0.4)
  expect_equal(pvalsmin_nosupport, mybtest[[1]]$p.value)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))

  pvals <- vapply(c("doubled", "kulinskaya", "minlikelihood"), function(m)
    twosidedpval(q = mybtest[[2]]$statistic, CDF = pbinom, continuous = FALSE,
                 method = m, Aloc = n[1] * 0.4, supportlim = c(0, n[1]),
                 size = n[1], prob = 0.4), NA_real_, USE.NAMES = FALSE)
  expect_equal(pvals[3], mybtest[[2]]$p.value)
  pvalsmin_nosupport <- twosidedpval(q = mybtest[[2]]$statistic, CDF = pbinom,
                 continuous = FALSE, method = "m", Aloc = n[1] * 0.4,
                 size = n[1], prob = 0.4)
  expect_equal(pvalsmin_nosupport, mybtest[[2]]$p.value)
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
})


# Poisson

test_that("twosidedpval works with Poisson distribution", {
  skip_on_cran()
  pvals <- vapply(c("doubled", "kulinskaya", "minlikelihood"), function(m)
    twosidedpval(q = 3, CDF = ppois, continuous = FALSE,
                 method = m, Aloc = 10, lambda = 10,
                 supportlim = c(0, Inf)), NA_real_, USE.NAMES = FALSE)
  pvalsmin_nosupport <- twosidedpval(q = 3, CDF = ppois, continuous = FALSE,
                 method = "m", Aloc = 10, lambda = 10, supportlim = c(0, Inf))
  lapply(pvals, function(p) expect_true(is.btwn01(p)))
  expect_equal(pvals[3], pvalsmin_nosupport)
})


