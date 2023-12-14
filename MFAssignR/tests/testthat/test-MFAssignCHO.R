patrick::with_parameters_test_that("MFAssignCHO works", {
#  load("../../data/Raw_Neg_ML.rda")
#   data <- system.file('data', 'Raw_Neg_ML.rda', package = "MFAssignR")
  load(system.file("data", "Raw_Neg_ML.rda",
                   package = "MFAssignR",
                   mustWork = TRUE
  ))
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

