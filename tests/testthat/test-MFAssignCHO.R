load_expected <- function(mode) {
  unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))
  ambig <- readRDS(file.path("test-data", paste0(mode, "_cho_ambig.rds")))
  none <- readRDS(file.path("test-data", paste0(mode, "_cho_none.rds")))
  return(list(Unambig = unambig, Ambig = ambig, None = none))
}

patrick::with_parameters_test_that("MFAssignCHO works",
  {
    peaks <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    actual <- MFAssignCHO(
      peaks = peaks$Mono,
      isopeaks = peaks$Iso,
      ionMode = mode
    )

    expected <- load_expected(mode)
    # keys <- c("Ambig", "Unambig", "None")
    # expect_equal(actual[keys], expected)

    expect_equal(actual$Ambig, expected$Ambig)
    expect_equal(actual$Unambig, expected$Unambig)
    expect_equal(actual$None, expected$None)
  },
  mode = c("pos", "neg")
)
