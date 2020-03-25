#context("multRM")

test_that("multRM with one dimension equals RM", {
  t1 <- multRM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", within = c("region", "feature"), iter = 1)
  t2 <- RM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", no.subf = 2, iter = 1)
  expect_equal(t1$WTS, t2$WTS)
})
