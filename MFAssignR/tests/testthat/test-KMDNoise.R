test_that("KMDNoise works", {
  load(system.file("data", "Raw_Neg_ML.rda",
    package = "MFAssignR",
    mustWork = TRUE
  ))
  names(Raw_Neg_ML) <- c("mass", "intensity")
  actual <- KMDNoise(Raw_Neg_ML)
  expect_equal(actual$Noise, 346.0706, tolerance = 1e-4)
})


test_that("SNplot works", {
  load(system.file("data", "Raw_Neg_ML.rda",
    package = "MFAssignR",
    mustWork = TRUE
  ))
  head(Raw_Neg_ML)
  Raw_Neg_ML <- dplyr::relocate(Raw_Neg_ML, intensity) |> dplyr::rename(mass = m.z)
  SNplot(Raw_Neg_ML, 1, 300, 100, 1000)
})