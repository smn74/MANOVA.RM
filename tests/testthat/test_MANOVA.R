#library(testthat)
library(MANOVA.RM)
context("MANOVA output")

test_that("MANOVA type",{
  data(EEG)
  mod1 <- MANOVA(resp ~ sex * diagnosis, data = EEG,
                 subject = "id", iter = 10, CPU = 1)
  expect_is(mod1, "MANOVA")
})


test_that("example 1: one-way",{
  mod <- MANOVA(resp ~ diagnosis, data = EEG, subject = "id", iter = 10, CPU = 1)
  names(mod$WTS) <- NULL
  expect_equal(mod$WTS[1], 53.553)
})

test_that("example 2: two-way",{
  mod <- MANOVA(resp ~ sex * diagnosis, data = EEG, subject = "id", iter = 10, CPU = 1)
  names(mod$WTS) <- NULL
  expect_equal(mod$WTS[1], 12.604)
})

test_that("example 3: three-way",{
  expect_warning(mod <- MANOVA(resp ~ sex * diagnosis * age, data = EEG, subject = "id", 
                               iter = 10, CPU = 1))
  names(mod$WTS) <- NULL
  expect_equal(mod$WTS[1], 28.874)
})

test_that("example 4: nested design",{
  library(GFD)
  data(curdies)
  set.seed(123)
  curdies$dug2 <- curdies$dugesia + rnorm(36)
  mod <- MANOVA.wide(cbind(dugesia, dug2) ~ season + season:site, data = curdies, iter = 10,
                       nested.levels.unique = TRUE, CPU = 1)
  names(mod$WTS) <- NULL
  expect_equal(mod$WTS[1], 6.999)
})


test_that("singular covariance matrix", {
  expect_warning(test <- MANOVA(resp ~ diagnosis*sex*age, data = EEG, iter = 10, subject = "id", CPU = 1))
})