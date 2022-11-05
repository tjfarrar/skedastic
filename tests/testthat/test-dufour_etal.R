context("dufour_etal works for two lm examples across all hypothesis tests in package")

htest <- c("bamset", "breusch_pagan", "carapeto_holt", "cook_weisberg",
           "diblasi_bowman", "evans_king", "glejser",
           "goldfeld_quandt", "harvey", "honda", "horn",
           "li_yao", "rackauskas_zuokas", "simonoff_tsai",
           "szroeter", "verbyla", "white", "zhou_etal")

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
  age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

# theargs <- formals(dufour_etal)


theargs <- lapply(1:length(htest), function(i) list("mainlm" = list(carslm, bostonlm)))

names(theargs) <- htest

# theargs$bickel$fitmethod <- c("lm", "rlm")
theargs$carapeto_holt$alternative <- c("greater", "less", "two.sided")
theargs$diblasi_bowman$distmethod <- c("moment.match", "bootstrap")
theargs$evans_king$method <- c("GLS", "LM")
theargs$goldfeld_quandt$alternative <- c("greater", "less", "two.sided")
theargs$honda$alternative <- c("greater", "less", "two.sided")
theargs$horn$restype <- c("ols", "blus")
theargs$horn$alternative <- c("greater", "less", "two.sided")
theargs$li_yao$method <- c("cvt", "alrt")
theargs$simonoff_tsai$method <- c("score")
theargs$zhou_etal$method <- c("pooled", "covariate-specific", "hybrid")

for (l in 1:length(theargs)) {
  theargs[[l]]$hettest <- names(theargs)[l]
}

allargs <- lapply(1:length(theargs), function(i) expand.grid(theargs[[i]],
                stringsAsFactors = FALSE))
names(allargs) <- names(theargs)

unlist(lapply(allargs, nrow))

allargs$zhou_etal <- allargs$zhou_etal[c(-4, -6), ]

# lapply(1:length(theargs), function(l)
#   test_that("linear regression works for dufour_etal with two regression models and each htest with method argument",
#   {pvals <- vapply(1:nrow(allargs[[l]]),
#     function(i) do.call(what = dufour_etal,
#              args = append(list("R" = 10L), unlist(allargs[[l]][i, ],
#                   recursive = FALSE)))$p.value, NA_real_)
#   lapply(pvals, function(p) expect_true(is.btwn01(p)))}))

lapply(1:length(theargs), function(l) {print(paste0("Test: ", names(allargs)[l]))
  test_that("dufour_etal with two regression models and each htest with method argument (normal errors)",
            {skip_on_cran()
              pvals <- vapply(1:nrow(allargs[[l]]),
            function(i) {print(paste0(i, " of ", nrow(allargs[[l]])))
            do.call(what = dufour_etal,
            args = append(list("R" = 10L), unlist(allargs[[l]][i, ],
            recursive = FALSE)))$p.value}, NA_real_)
            lapply(pvals, function(p) expect_true(is.btwn01(p)))}) } )


lapply(1:length(theargs), function(l) {print(paste0("Test: ", names(allargs)[l]))

  test_that("dufour_etal with two regression models and each htest with method argument (t errors)",
            {skip_on_cran()
            pvals <- vapply(1:nrow(allargs[[l]]),
            function(i) {print(paste0(i, " of ", nrow(allargs[[l]])))
            do.call(what = dufour_etal,
            args = append(list("R" = 10L, "errorgen" = stats::rt,
                    "errorparam" = list("df" = 3)), unlist(allargs[[l]][i, ],
            recursive = FALSE)))$p.value}, NA_real_)
            lapply(pvals, function(p) expect_true(is.btwn01(p)))}) } )


lapply(1:length(theargs), function(l) {print(paste0("Test: ", names(allargs)[l]))
  test_that("dufour_etal with two regression models and each htest with method argument (uniform errors)",
            {skip_on_cran()
            pvals <- vapply(1:nrow(allargs[[l]]),
            function(i) {print(paste0(i, " of ", nrow(allargs[[l]])))
            do.call(what = dufour_etal,
            args = append(list("R" = 10L, "errorgen" = stats::runif,
                    "errorparam" = list("min" = -1, "max" = 1)), unlist(allargs[[l]][i, ],
            recursive = FALSE)))$p.value}, NA_real_)
            lapply(pvals, function(p) expect_true(is.btwn01(p)))}) } )
