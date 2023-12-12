test_that("Isotope filtering works", {
  load(system.file("data", "Raw_Neg_ML.rda",
    package = "MFAssignR",
    mustWork = TRUE
  ))
  expected <- readRDS("test-data/isotopes.rda")

  actual <- IsoFiltR(Raw_Neg_ML)

  expect_equal(actual, expected)
})
