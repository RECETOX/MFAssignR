test_that("Recal works when isopeaks is none", {
    unambig <- readRDS(file.path("test-data/pos_cho_unambig.rds"))
    mono <- readRDS(file.path("test-data/pos_recal_mono.rds"))
    mono <- mono[, c(1:2)]

  expected <- readRDS(file.path("test-data", "recal_2cols_isopeaks-none.rds"))
  actual <- MFAssignR::Recal(
      unambig,
      peaks = mono,
      isopeaks = "none",
      mzRange = 80,
      mode = "pos",
      SN = 0,
      series1 = "O4_H_11",
      series2 = "O2_H_6",
      series3 = "O7_H_13",
      series4 = "O_H_13",
      series5 = "O4_H_6",
      series6 = NA,
      series7 = NA,
      series8 = NA,
      series9 = NA,
      series10 = NA
    )

  saveRDS(actual, "test-data/recal-actual.rds")

  expect_equal(actual$Mono, expected$Mono)
  expect_equal(actual$Iso, expected$Iso)
  expect_equal(actual$RecalList, expected$RecalList)
})