test_that("MFAssignCHO works", {
  data <- readRDS("test-data/isotopes.rda")
  expected <- readRDS("test-data/mfassignCHO.rda")
  
  actual <- MFAssignCHO(peaks = data$Mono, isopeaks = data$Iso, ionMode = "neg")
  
  expect_equal(actual$Unambig, expected$Unambig)
  expect_equal(actual$Ambig, expected$Ambig)
  expect_equal(actual$None, expected$None)
})

#saveRDS(actual, "test-data/mfassignCHO.rda")


# unambig <- expected$Unambig
# ambig <- expected$Ambig
# none <- expected$none
# 
# saveRDS(unambig, "test-data/unambig.rda")
# saveRDS(ambig, "test-data/ambig.rda")
# saveRDS(none, "test-data/none.rda")
# 
