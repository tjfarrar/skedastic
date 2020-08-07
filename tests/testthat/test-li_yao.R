context("li_yao works for two lm examples across all argument permutations")

# theargs <- formals(li_yao)

test_that("linear regression works with all combinations of formals", {
  carslm <- lm(dist ~ speed, data = cars)
  bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
                   age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)
  theargs <- list("method" = c("cvt", "alrt"), baipanyin = c(TRUE, FALSE),
              "mainlm" = list(carslm, bostonlm))
  allargs <- expand.grid(theargs, stringsAsFactors = FALSE)
  allargs <- allargs[-which(allargs$method == "alrt" & allargs$baipanyin), ]

  pvals <- vapply(1:nrow(allargs), function(i) do.call(what = li_yao,
              args = append(list("statonly" = FALSE),
              unlist(allargs[i, ], recursive = FALSE)))$p.value, NA_real_)
  lapply(1:length(pvals), function(i) expect_true(is.btwn01(pvals[i])))
})
