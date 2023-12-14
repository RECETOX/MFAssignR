load_test_data <- function() {
  iso <- readRDS("test-data/iso_neg.rda")
  mono <- readRDS("test-data/mono_neg.rda")
  mfassingr <- readRDS("test-data/mfassign_Unambig_neg.rda")

  list(iso = iso, mono = mono, mfassingr = mfassingr)
}

patrick::with_parameters_test_that("Recal works", {
  data <- load_test_data()

  actual <- MFAssignR::Recal(
    data$mfassingr,
    peaks = data$mono,
    isopeaks = data$iso,
    mzRange = 30,
    mode = mode,
    SN = 0,
    series1 = "O10_H_10",
    series2 = "O5_H_6",
    series3 = "O7_H_8",
    series4 = "O8_H_8",
    series5 = "O4_H_2",
    series6 = "O10_H_9"
  )
  recal <- readRDS("test-data/recal_neg.rda")
  expect_equal(actual, recal)
},
  mode = c("neg"),
)
test_that("RecalList works", {
  data <- load_test_data()
  recalist <- readRDS("test-data/recalist_neg.rda")

  actual <- MFAssignR::RecalList(data$mfassingr)

  expect_equal(actual, recalist, tolerance = 1e-2)
})
