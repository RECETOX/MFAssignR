load_test_data <- function() {
  iso <- readRDS("test-data/iso.rda")
  mono <- readRDS("test-data/mono.rda")
  mfassingr <- readRDS("test-data/mfassignR.rda")

  list(iso = iso, mono = mono, mfassingr = mfassingr)
}

patrick::with_parameters_test_that("Recal works", {
  data <- load_test_data()

  actual <- MFAssignR::Recal(
    data$mfassingr[["Unambig"]],
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
  recal <- readRDS(paste0("test-data/recal-", mode, ".rda"))
  expect_equal(actual, recal)
},
  mode = c("neg", "pos"),
)

test_that("RecalList works", {
  data <- load_test_data()
  recalist <- readRDS("test-data/recalist.rda")

  actual <- MFAssignR::RecalList(data$mfassingr[["Unambig"]])

  expect_equal(actual, recalist, tolerance = 1e-2)
})
