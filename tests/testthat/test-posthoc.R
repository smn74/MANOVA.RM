context("post-hoc")

test_that("post-hoc without interaction equals simple model", {
  t1 <- MANOVA.wide(cbind(brainrate_central, brainrate_frontal) ~ sex*diagnosis, data = EEGwide, iter = 10, seed = 123)
  t2 <- MANOVA.wide(cbind(brainrate_central, brainrate_frontal) ~ diagnosis, data = EEGwide, iter = 10, seed = 123)
  expect_equal(simCI(t1, contrast = "pairwise", type = "Tukey", interaction = FALSE, factor = "diagnosis"), 
  simCI(t2, contrast = "pairwise", type = "Tukey"))
})
