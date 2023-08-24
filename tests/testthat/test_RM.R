library(testthat)
library(MANOVA.RM)
context("RM output")

test_that("Means correctly calculated",{
  if(requireNamespace("data.table")){
    library(data.table)
  data("o2cons")
  ox <- as.data.table(o2cons)
  m <- ox[, mean(O2), by = .(Group, Time, Staphylococci)]
  t1 <- RM(O2 ~ Group * Time * Staphylococci, data = o2cons, subject = "Subject",
           iter = 1, within = c("Time", "Staphylococci"))
  t2 <- RM(O2 ~ Group * Staphylococci * Time , data = o2cons, subject = "Subject",
           iter = 1, within = c("Time", "Staphylococci"))
  expect_equal(round(m[Group == "P" & Time == 6 & Staphylococci == 0, V1], 3),
               t1$Descriptive[1, "Means"], t2$Descriptive[1, "Means"])
  expect_equal(round(m[Group == "P" & Time == 18 & Staphylococci == 0, V1], 3),
               t1$Descriptive[5, "Means"], t2$Descriptive[3, "Means"])
}
  })

test_that("example 1: 1 whole, 2 sub",{
  oxy <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
            subject = "Subject", no.subf = 2, iter = 100, CPU = 1)
  expect_equal(oxy$WTS[1], 11.167)
})

test_that("example 2: 2 whole, 2 sub", {
  EEG_model <- RM(resp ~ sex * diagnosis * feature * region,
                  data = EEG, subject = "id", no.subf = 2, resampling = "WildBS",
                  iter = 100, CPU = 1)
  expect_equal(EEG_model$WTS[1], 9.973)
})

test_that("example 3: only subplot factors",{
  data("o2cons")
  oxy <- o2cons[o2cons$Group == "P", ]
  submodel <-RM(O2~ Staphylococci * Time, data = oxy, 
                within = c("Staphylococci", "Time"), 
                iter = 2, CPU = 1, subject = "Subject")
  expect_equal(submodel$WTS[1], 2.69)
})

test_that("correct mean calculation",{
  if(requireNamespace("GFD")){
    library(GFD)
  data("o2cons")
  oxy <- o2cons[o2cons$Group == "P" & o2cons$Staphylococci == 1, ]
  mod1 <- RM(O2 ~ Time, data = oxy, subject = "Subject", no.subf = 1,
             iter = 10, CPU = 1)
  mod2 <- GFD(O2 ~ Time, data = oxy, nperm = 10)
  # WTS and ATS differ, since GFD ignores dependency structure
  expect_equal(c(mod1$Descriptive["Means"]), c(round(mod2$Descriptive["Means"], 3)))
}
  })

test_that("wrong interaction",{
  data(EEG)
  expect_error(RM(resp ~ sex + feature * region, data = EEG,
                  subject = "id", iter = 10, CPU = 1, no.subf = 2))
})

test_that("missing values", {
  data(EEG)
  EEG2 <- EEG[-5, ]
  expect_error(RM(resp ~ sex * feature * region, data = EEG2, no.subf = 2,
                  subject = "id", iter = 1, CPU = 1))
})



test_that("error if within not last in formula", {
  data(EEG)
  expect_error(RM(O2 ~ Staphylococci * Time * Group, data = o2cons, subject = "Subject",
                  iter = 1, within = c("Time", "Staphylococci")))
})

