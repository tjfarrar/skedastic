context("pRQF functions work for various choices of A, B, Sigma")

# Should error

n1 <- 10
A1 <- matrix(data = 1, nrow = n1, ncol = n1)
B1 <- matrix(data = 0, nrow = n1, ncol = n1)
Bneg <- matrix(data = -1, nrow = n1, ncol = n1)

set.seed(12345)
Arnd0 <- matrix(data = runif(n1 ^ 2, min = -1, max = 1),
                nrow = n1, ncol = n1)
Arnd <- matrix(data = runif(n1 ^ 2), nrow = n1, ncol = n1)
Arnd <- Arnd %*% t(Arnd)
Brnd <- matrix(data = runif(n1 ^ 2), nrow = n1, ncol = n1)
Brnd <- Brnd %*% t(Brnd)

test_that("pRQF: test some erroneous settings", {

  expect_error(pRQF(r = 1, A = A1, B = B1))
  expect_error(pRQF(r = 1, A = A1, B = Bneg))
  expect_error(pRQF(r = 1, A = Arnd0, B = Arnd0))
})

theeigen <- eigen(solve(Brnd) %*% Arnd, only.values = TRUE)$values
minsupport <- min(theeigen)
maxsupport <- max(theeigen)
test_that("pRQF: test some trivial cases", {
  expect_equal(pRQF(r = 1, A = Arnd, B = Arnd), 0.5)
  expect_equal(pRQF(r = 0, A = Arnd, B = Arnd), 0)
  expect_equal(pRQF(r = 1, A = B1, B = A1), 1)
  expect_equal(pRQF(r = minsupport, A = Arnd, B = Brnd), 0)
  expect_equal(suppressWarnings(pRQF(r = maxsupport, A = Arnd, B = Brnd)), 1)
  expect_equal(pRQF(r = 1.5, A = diag(n1), B = A1, algorithm = "imhof"),
               pRQF(r = 1.5, A = diag(n1), B = A1, algorithm = "integrate"))
})

set.seed(12345)
Abinrnd <- matrix(data = sample(c(0, -1, 1), n1 ^ 2, replace = TRUE),
                  nrow = n1, ncol = n1)
Abinrnd <- Abinrnd %*% t(Abinrnd)
Bbinrnd <- matrix(data = sample(c(0, -1, 1), n1 ^ 2, replace = TRUE),
                  nrow = n1, ncol = n1)
Bbinrnd <- Bbinrnd %*% t(Bbinrnd)

theeigen2 <- eigen(solve(Bbinrnd) %*% Abinrnd, only.values = TRUE)$values
minsupport2 <- min(theeigen2)
maxsupport2 <- max(theeigen2)
mySigma <- Arnd0 %*% t(Arnd0)

test_that("pRQF: test some trivial cases", {
  expect_equal(pRQF(r = 1, A = Abinrnd, B = Abinrnd), 0.5)
  expect_equal(pRQF(r = 0, A = Abinrnd, B = Abinrnd), 0)
  expect_equal(pRQF(r = minsupport2, A = Abinrnd, B = Bbinrnd), 0)
  expect_equal(suppressWarnings(pRQF(r = maxsupport2, A = Abinrnd, B = Bbinrnd)), 1)
  expect_equal(pRQF(r = exp(1), A = Arnd, B = Brnd, Sigma = mySigma,
                    algorithm = "integrate"),
               pRQF(r = exp(1), A = Arnd, B = Brnd, Sigma = mySigma,
                    algorithm = "imhof"))
})


test_that("pRQF: Monte Carlo results", {
  skip_on_cran()
  ### Monte Carlo
  S <- 1e4
  set.seed(1234)
  eps <- replicate(S, rnorm(n1), simplify = FALSE)
  eps2 <- replicate(S, mvtnorm::rmvnorm(n = 1, sigma = mySigma), simplify = FALSE)
  R <- unlist(lapply(eps, function(x)
    (t(x) %*% Arnd %*% x) / (t(x) %*% Brnd %*% x)))
  R2 <- unlist(lapply(eps2, function(x)
    (x %*% Arnd %*% t(x)) / (x %*% Brnd %*% t(x))))

  expect_true(abs(sum(R < 1) / S - pRQF(r = 1, A = Arnd, B = Brnd)) <= 1e-2)
  expect_true(abs(sum(R2 < 0.875) / S -
              pRQF(r = 0.875, A = Arnd, B = Brnd, Sigma = mySigma)) <= 1e-2)
  expect_true(abs(sum(R >= 1.1) / S - pRQF(r = 1.1, A = Arnd, B = Brnd,
              lower.tail = FALSE)) <= 1e-2)
  expect_true(abs(sum(R2 >= 0.62) / S -
              pRQF(r = 0.62, A = Arnd, B = Brnd, Sigma = mySigma,
                   lower.tail = FALSE)) <= 1e-2)
})
