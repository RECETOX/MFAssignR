expected_paths <- function(mode) {
  return(list(
    unambig  = file.path("test-data", paste0(mode, "_unambig.rds")),
    ambig  = file.path("test-data", paste0(mode, "_ambig.rds")),
    none  = file.path("test-data", paste0(mode, "_none.rds"))
  ))
}

load_expected <- function(mode) {
  paths <- expected_paths(mode)
  return(list(
    Unambig = readRDS(paths$unambig),
    Ambig = readRDS(paths$ambig),
    None = readRDS(paths$none)
  ))
}

save_expected <- function(actual, mode) {
  paths <- expected_paths(mode)
  saveRDS(actual$Unambig, paths$unambig)
  saveRDS(actual$Ambig, paths$ambig)
  saveRDS(actual$None, paths$none)
}

patrick::with_parameters_test_that("MFAssign works",
  {
    peaks <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    actual <- MFAssign(
      peaks = peaks$Mono,
      isopeaks = peaks$Iso,
      ionMode = mode
    )
    save_expected(actual, mode)
    expected <- load_expected(mode)

    expect_equal(actual$Ambig, expected$Ambig)
    expect_equal(actual$Unambig, expected$Unambig)
    expect_equal(actual$None, expected$None)
  },
  mode = c("pos", "neg")
)
