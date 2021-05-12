library(MANOVA.RM)
context("MANOVA vs. MANOVA.wide")

test_that("MANOVA equals MANOVA.wide",{
  data(EEG)
  data("EEGwide")
  mod1 <- MANOVA(resp ~ sex * diagnosis, data = EEG,
                 subject = "id", iter = 10, CPU = 1)
  mod2 <- MANOVA.wide(cbind(brainrate_temporal, brainrate_frontal, brainrate_central, complexity_temporal, 
                            complexity_frontal, complexity_central) ~ sex * diagnosis,
                          data = EEGwide, iter = 10, CPU = 1)
  expect_equal(mod1$WTS, mod2$WTS)
})


test_that("one dimension", {
  object <- MANOVA.wide(brainrate_central ~ diagnosis, data = EEGwide,
                        iter = 1000, CPU = 1)
  expect_equal(object$MATS[1], 41.75)
})