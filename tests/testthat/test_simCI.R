library(testthat)
library(MANOVA.RM)
context("Post-hoc tests")

test_that("example 1: user-defined", {
object <- MANOVA.wide(cbind(brainrate_temporal, brainrate_central) ~ diagnosis, data = EEGwide,
                      iter = 1000, CPU = 1)
H <- as.matrix(cbind(rep(1, 5), -1*Matrix::Diagonal(5)))
pp <- simCI(object, contrast = "user-defined", contmat = H, silent = TRUE)
expect_equal(pp$Estimate[1], 0.013)
})

test_that("example 2: Dunnett", {
  test <- MANOVA(resp ~ diagnosis, data = EEG, iter = 10, subject = "id", CPU = 1)
  pp <- simCI(test, contrast = "pairwise", type = "Dunnett", silent = TRUE)
  expect_equal(pp$Estimate[1], 1.798)
})

test_that("example 3: 3-way", {
  expect_warning(test <- MANOVA(resp ~ diagnosis*sex*age, data = EEG, iter = 10, 
                                subject = "id", CPU = 1))
  pp <- simCI(test, contrast = "pairwise", type = "Dunnett", silent = TRUE)
  expect_equal(pp$Estimate[1], 4.196)
})

test_that("nested design", {
  if(requireNamespace("GFD")){
    library(GFD)
  data(curdies)
  set.seed(123)
  curdies$dug2 <- curdies$dugesia + rnorm(36)
  fit1 <- MANOVA.wide(cbind(dugesia, dug2) ~ season + season:site, data = curdies, iter = 10,
                      nested.levels.unique = TRUE, seed = 123, CPU = 1)
  expect_error(simCI(fit1, contrast = "pairwise", type = "Dunnett"))
}
  })