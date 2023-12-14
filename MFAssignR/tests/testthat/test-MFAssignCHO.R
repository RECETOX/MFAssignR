patrick::with_parameters_test_that("MFAssignCHO works", {
  load(system.file("data", "Raw_Neg_ML.rda",
                   package = "MFAssignR",
                   mustWork = TRUE
  ))
  data <- IsoFiltR(Raw_Neg_ML)
  actual <- MFAssignCHO(
    peaks = data$Mono,
    isopeaks = data$Iso,
    ionMode = mode
  )
  expected <- readRDS(paste0("test-data/mfassignCHO_testData-", mode, ".rda"))
  expect_equal(actual$Ambig, expected$Ambig)
  expect_equal(actual$Unambig, expected$Unambig)
  expect_equal(actual$None, expected$None)
},
mode = c("neg", "pos")
)
