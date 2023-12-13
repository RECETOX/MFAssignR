patrick::with_parameters_test_that("MFAssignCHO works", {
  load("../../inst/data/Raw_Neg_ML.rda")
  data <- MFAssignR::IsoFiltR(Raw_Neg_ML)
  actual <- MFAssignR::MFAssignCHO(
    peaks = data$Mono,
    isopeaks = data$Iso,
    ionMode = mode
  )
  expected <- readRDS(paste0("test-data/mfassignCHO_testData-", mode, ".rda"))
  expect_equal(actual, expected)
},
mode = c("neg", "pos")
)

