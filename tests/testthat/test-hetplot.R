context("hetplot works for two lm examples across all argument permutations")

# theargs <- formals(rackauskas_zuokas)

carslm <- lm(dist ~ speed, data = cars)
bostonlm <- lm(medv ~ crim + zn + indus + chas + nox + rm +
            age + dis + rad + tax + ptratio + b + lstat, data = BostonHousing)

theargs <- vector("list", 2)

theargs <- lapply(list(carslm, bostonlm), function(mod)
                lapply(c(NA, "png", "bmp", "jpeg", "tiff"),
                function(ft) list("horzvar" = c("index", "explanatory", "fitted.values",
                "fitted.values2", "log_explanatory", "speed", "log_speed"),
                "vertvar" = c("res", "res_blus",
                "res_stand", "res_constvar", "res_stud"), "vertfun" = c("2",
                "3", "abs", "identity", "sin"), "filetype" = ft, "values" = c(TRUE, FALSE),
                "mainlm" = list(mod))))

theargs <- unlist(theargs, recursive = FALSE)

allargs <- lapply(theargs, function(ta) expand.grid(ta, stringsAsFactors = FALSE))

for (i in 6:10) {
  allargs[[i]]$horzvar[allargs[[i]]$horzvar == "speed"] <- "crim"
  allargs[[i]]$horzvar[allargs[[i]]$horzvar == "log_speed"] <- "log_crim"
}

# filetypes <- vapply(allargs, function(a) a$filetype[1], NA_character_)

# do.call(what = hetplot,
#         args = append(list(), unlist(allargs[[1]][1, ], recursive = FALSE)))
#

## par(new) causing some problems

# hetplot(mainlm = carslm, horzvar = "index", vertvar = "res", vertfun = "sin",
#  filetype = NA, values = TRUE)
#
# hetplot(mainlm = carslm, horzvar = "index", vertvar = "res", vertfun = "abs",
#         filetype = NA, values = TRUE)
#

test_that("hetplot works with all combinations of formals for carslm, filetype NA", {
  doplots <- lapply(1:nrow(allargs[[1]]), function(i) do.call(what = hetplot,
  args = append(list(), unlist(allargs[[1]][i, ], recursive = FALSE))))
  expect_true(length(doplots) == 350L)
})

test_that("hetplot works with all combinations of formals for carslm, filetype png", {
  doplots <- lapply(1:nrow(allargs[[2]]), function(i) do.call(what = hetplot,
            args = append(list(), unlist(allargs[[2]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "png")
  expect_true(length(lpng) == sum(allargs[[2]]$filetype == "png", na.rm = TRUE))
  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for carslm, filetype bmp", {
  skip_on_cran()
  doplots <- lapply(1:nrow(allargs[[3]]), function(i) do.call(what = hetplot,
            args = append(list(), unlist(allargs[[3]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "bmp")
  expect_true(length(lpng) == sum(allargs[[3]]$filetype == "bmp", na.rm = TRUE))
  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for carslm, filetype jpeg", {
  skip_on_cran()
  doplots <- lapply(1:nrow(allargs[[4]]), function(i) do.call(what = hetplot,
            args = append(list(), unlist(allargs[[4]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "jpeg")
  expect_true(length(lpng) == sum(allargs[[4]]$filetype == "jpeg", na.rm = TRUE))
  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for carslm, filetype tiff", {
  skip_on_cran()
  doplots <- lapply(1:nrow(allargs[[5]]), function(i) do.call(what = hetplot,
              args = append(list(), unlist(allargs[[5]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "tiff")
  expect_true(length(lpng) == sum(allargs[[5]]$filetype == "tiff", na.rm = TRUE))
  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

allargs[[6]] <- allargs[[6]][-which(allargs[[6]]$horzvar %in%
                            c("explanatory", "log_explanatory")), ]
test_that("hetplot works with all combinations of formals for bostonlm, filetype NA", {
  skip_on_cran()
  doplots <- lapply(1:nrow(allargs[[6]]), function(i) do.call(what = hetplot,
      args = append(list(), unlist(allargs[[6]][i, ], recursive = FALSE))))
  expect_true(length(doplots) == 250L)
})

test_that("hetplot works with all combinations of formals for bostonlm, filetype png", {
  skip_on_cran()

  doplots <- lapply(1:nrow(allargs[[7]]), function(i) do.call(what = hetplot,
      args = append(list(), unlist(allargs[[7]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "png")

  numall <- sum((allargs[[7]]$horzvar == "explanatory" |
                  allargs[[7]]$horzvar == "log_explanatory") &
                  allargs[[7]]$filetype == "png", na.rm = TRUE)
  expectedfiles <- numall * 12 + sum(allargs[[7]]$filetype == "png", na.rm = TRUE)
  expect_true(length(lpng) == expectedfiles)
  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for bostonlm, filetype bmp", {
  skip_on_cran()

  doplots <- lapply(1:nrow(allargs[[8]]), function(i) do.call(what = hetplot,
      args = append(list(), unlist(allargs[[8]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "bmp")

  numall <- sum((allargs[[8]]$horzvar == "explanatory" |
                   allargs[[8]]$horzvar == "log_explanatory") &
                  allargs[[8]]$filetype == "bmp", na.rm = TRUE)
  expectedfiles <- numall * 12 + sum(allargs[[8]]$filetype == "bmp", na.rm = TRUE)
  expect_true(length(lpng) == expectedfiles)

  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for bostonlm, filetype jpeg", {
  skip_on_cran()

  doplots <- lapply(1:nrow(allargs[[9]]), function(i) do.call(what = hetplot,
      args = append(list(), unlist(allargs[[9]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "jpeg")

  numall <- sum((allargs[[9]]$horzvar == "explanatory" |
                   allargs[[9]]$horzvar == "log_explanatory") &
                  allargs[[9]]$filetype == "jpeg", na.rm = TRUE)
  expectedfiles <- numall * 12 + sum(allargs[[9]]$filetype == "jpeg", na.rm = TRUE)
  expect_true(length(lpng) == expectedfiles)

  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})

test_that("hetplot works with all combinations of formals for bostonlm, filetype tiff", {
  skip_on_cran()

  doplots <- lapply(1:nrow(allargs[[10]]), function(i) do.call(what = hetplot,
      args = append(list(), unlist(allargs[[10]][i, ], recursive = FALSE))))
  lpng <- list.files(paste0(tempdir(), "/hetplot/"), pattern = "tiff")

  numall <- sum((allargs[[10]]$horzvar == "explanatory" |
                   allargs[[10]]$horzvar == "log_explanatory") &
                  allargs[[10]]$filetype == "tiff", na.rm = TRUE)
  expectedfiles <- numall * 12 + sum(allargs[[10]]$filetype == "tiff", na.rm = TRUE)
  expect_true(length(lpng) == expectedfiles)

  lapply(paste0(tempdir(), "/hetplot/", lpng), unlink)
})
