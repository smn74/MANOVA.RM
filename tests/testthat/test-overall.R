#context("Overall functioning")

test_that("RM-function", {
  oxy <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons,
            subject = "Subject", no.subf = 2, iter = 1000, resampling = "Perm", seed = 123)
  expect_equal_to_reference(summary(oxy), file = "RM.rds")
})


test_that("MANOVA",{
  EEG_mod <- MANOVA(resp ~ sex * diagnosis, data = EEG, subject = "id", resampling = "paramBS", seed = 123)
  expect_equal_to_reference(summary(EEG_mod), file = "MANOVA.rds")
})

test_that("MANOVAwide", {
  t1 <- MANOVA.wide(cbind(brainrate_central, brainrate_frontal) ~ sex*diagnosis, data = EEGwide, seed = 123)
  expect_equal_to_reference(summary(t1), file = "MANOVAwide.rds")
})

test_that("multRM", {
  library(tidyr)
  eeg <- spread(EEG, feature, resp)
  fit <- multRM(cbind(brainrate, complexity) ~ sex + region, data = eeg, subject = "id", within = "region", seed = 123)
  expect_equal_to_reference(summary(fit), file = "multRM.rds")
})

