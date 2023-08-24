library(testthat)
library(MANOVA.RM)
context("multRM")

test_that("Means correctly calculated",{
  if(requireNamespace("data.table")){
    library(data.table)
  data("EEG")
  eeg2 <- as.data.table(EEG)
  m <- eeg2[, mean(resp), by = .(diagnosis, region, feature)]
  t1 <- multRM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", 
               within = c("region", "feature"), iter = 1, CPU =1)
  t2 <- multRM(resp ~ diagnosis * feature * region, data = EEG, subject = "id", 
               within = c("region", "feature"), iter = 1, CPU =1)
  expect_equal(round(m[diagnosis == "AD" & region == "central" & feature == "brainrate", V1], 3),
               t1$Descriptive[1, "resp"], t2$Descriptive[1, "resp"])
  expect_equal(round(m[diagnosis == "SCC" & region == "frontal" & feature == "brainrate", V1], 3),
               t1$Descriptive[15, "resp"], t2$Descriptive[14, "resp"])
}
  })


test_that("multRM equals RM",{
  data(EEG)
  t1 <- multRM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", 
               within = c("region", "feature"), iter = 1, CPU =1)
  t2 <- RM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", 
               no.subf = 2, iter = 1, CPU = 1)
  expect_equal(t1$WTS, t2$WTS)
})


test_that("missing values", {
  data(EEG)
  EEG2 <- EEG[-5, ]
  expect_error(multRM(resp ~ sex * feature * region, data = EEG2,
                      within = c("feature", "region"),
                  subject = "id", iter = 1, CPU = 1))
})


test_that("error if within not last in formula", {
  data(EEG)
  expect_error(multRM(resp ~ feature * region * diagnosis, data = EEG, subject = "id", 
                      within = c("region", "feature"), iter = 1, CPU =1))
})

