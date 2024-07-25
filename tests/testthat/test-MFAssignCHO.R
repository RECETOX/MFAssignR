make_path <- function(mode, identifier) {
  file.path("test-data", paste0(mode, "_cho_", identifier, ".rds"))
}

load_expected <- function(mode) {
  unambig <- readRDS(make_path(mode, "unambig"))
  ambig <- readRDS(make_path(mode, "ambig"))
  none <- readRDS(make_path(mode, "none"))
  return(list(Unambig = unambig, Ambig = ambig, None = none))
}

update_expected <- function(actual, mode) {
  saveRDS(actual$Ambig, file = make_path(mode, "ambig"))
  saveRDS(actual$Unambig, file = make_path(mode, "unambig"))
  saveRDS(actual$None, file = make_path(mode, "none"))
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

    expect_equal(actual$Ambig, expected$Ambig)
    expect_equal(actual$Unambig, expected$Unambig)
    expect_equal(actual$None, expected$None)
  },
  mode = c("pos", "neg")
)
