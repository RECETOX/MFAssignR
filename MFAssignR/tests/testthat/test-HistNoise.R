
test_that("HistNoise works", {
   load(system.file("data", "Raw_Neg_ML.rda",
    package = "MFAssignR",
    mustWork = TRUE
  ))
   actual <- HistNoise(Raw_Neg_ML)
   expect_equal(actual$Noise, 317.3483, tolerance = 1e-4)
 })
