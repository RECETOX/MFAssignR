test_that("HistNoise works", {
  load(file.path("test-data", "Raw_Neg_ML.rda"))
  actual <- HistNoise(Raw_Neg_ML)
  expect_equal(actual$Noise, 317.3483, tolerance = 1e-4)
})
