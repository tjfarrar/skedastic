context("alvm.fit works for cars lm and boston lm across several argument permutations")

set.seed(234)
carslm <- lm(dist ~ speed, data = cars)

myalvms <- list(alvm.fit(carslm, model = "linear", varselect = "qgcv.linear"),
                alvm.fit(carslm, model = "cluster", nclust = "elbow.mwd",
                         varselect = "hettest", hettest = "breusch_pagan"),
                alvm.fit(carslm, model = "polynomial"))

myanlvms <- list(anlvm.fit(carslm, g = function(x) x ^ 2,
                           varselect = "qgcv.linear"),
                anlvm.fit(carslm, g = function(x) x ^ 2, cluster = TRUE,
                         varselect = "hettest", hettest = "breusch_pagan"))


theargs1 <- list("bootsampmethod" = c("pairs", "wild"),
                 "bootCImethod" = c("pct", "bca", "stdnorm"),
                 "expand" = c(TRUE, FALSE),
                 rm_on_constraint = TRUE,
                 jackknife_point = FALSE,
                 Brequired = 5L, Bextra = 5L)

theargs2 <- list("bootsampmethod" = c("pairs", "wild"),
                 "bootCImethod" = c("pct", "bca", "stdnorm"),
                 "expand" = c(TRUE, FALSE),
                 rm_nonconverged = FALSE,
                 jackknife_point = FALSE,
                 Brequired = 5L, Bextra = 1L)


allargs1 <- expand.grid(theargs1, stringsAsFactors = FALSE)
allargs2 <- expand.grid(theargs2, stringsAsFactors = FALSE)

test_that("carslm ALVMs work with all combinations of formals", {
  skip_on_cran()
 alvm.cis <- lapply(1:length(myalvms), function(j) {
   # print(paste0("j = ", j))
   lapply(1:nrow(allargs1), function(i) {
     # print(paste0("i = ", i))
     do.call(what = avm.ci,
             args = append(list("object" = myalvms[[j]]), allargs1[i, ]))
   })
 })
 lapply(1:length(myalvms), function(j) {
   lapply(1:nrow(allargs1), function(i) {
     expect_true(all(is.finite(alvm.cis[[j]][[i]]$climits)))
   })
 })

})

test_that("carslm ANLVMs work with all combinations of formals", {
  skip_on_cran()
  anlvm.cis <- lapply(1:length(myanlvms), function(j) {
   # print(paste0("j = ", j))
   lapply(1:nrow(allargs2), function(i) {
     # print(paste0("i = ", i))
     do.call(what = avm.ci,
             args = append(list("object" = myanlvms[[j]]), allargs2[i, ]))
   })
  })
  lapply(1:length(myanlvms), function(j) {
   lapply(1:nrow(allargs2), function(j) {
     expect_true(all(is.finite(anlvm.cis[[j]][[i]]$climits)))
   })
  })

})

theargs3 <- list(rm_on_constraint = FALSE,
                 Brequired = 30L, Bextra = 30L,
                 retune = c(TRUE, FALSE))

theargs4 <- list(rm_nonconverged = FALSE,
                 Brequired = 30L, Bextra = 30L,
                 retune = c(TRUE, FALSE))

allargs3 <- expand.grid(theargs3, stringsAsFactors = FALSE)
allargs4 <- expand.grid(theargs4, stringsAsFactors = FALSE)


test_that("carslm ALVMs and ANLVMs work with and without retune", {
  skip_on_cran()

  alvm.cis <- lapply(1:length(myalvms), function(j) {
    # print(paste0("j = ", j))
    lapply(1:nrow(allargs3), function(i) {
      # print(paste0("i = ", i))
      do.call(what = avm.ci,
              args = append(list("object" = myalvms[[j]]), allargs3[i, ]))
    })
  })
  lapply(1:length(myalvms), function(j) {
    lapply(1:nrow(allargs3), function(i) {
      expect_true(all(is.finite(alvm.cis[[j]][[i]]$climits)))
    })
  })

  anlvm.cis <- lapply(1:length(myanlvms), function(j) {
    # print(paste0("j = ", j))
    lapply(1:nrow(allargs4), function(i) {
      # print(paste0("i = ", i))
      do.call(what = avm.ci,
              args = append(list("object" = myanlvms[[j]]), allargs4[i, ]))
    })
  })
  lapply(1:length(myanlvms), function(j) {
    lapply(1:nrow(allargs2), function(j) {
      expect_true(all(is.finite(anlvm.cis[[j]][[i]]$climits)))
    })
  })

})
