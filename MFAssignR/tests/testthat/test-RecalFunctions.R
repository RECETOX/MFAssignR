load_test_data <- function() {
  iso <- readRDS("test-data/iso.rda")
  mono <- readRDS("test-data/mono.rda")
  mfassingr <- readRDS("test-data/mfassignR.rda")

  list(iso = iso, mono = mono, mfassingr = mfassingr)
}

test_that("Recal works", {
  data <- load_test_data()
  recal <- readRDS("test-data/recal.rda")

  actual <- MFAssignR::Recal(data$mfassingr[["Unambig"]],
    peaks = data$mono,
    isopeaks = data$iso, mzRange = 30, mode = "neg", SN = 0,
    series1 = "O10_H_10", series2 = "O5_H_6", series3 = "O7_H_8",
    series4 = "O8_H_8", series5 = "O4_H_2", series6 = "O10_H_9"
  )

  expect_equal(actual$Mono, recal$Mono, tolerance = 1e-2)
  expect_equal(actual$Iso, recal$Iso, tolerance = 1e-2)
})

test_that("RecalList works", {
  data <- load_test_data()
  recalist <- readRDS("test-data/recalist.rda")

  actual <- MFAssignR::RecalList(data$mfassingr[["Unambig"]])

  expect_equal(actual, recalist, tolerance = 1e-2)
})
