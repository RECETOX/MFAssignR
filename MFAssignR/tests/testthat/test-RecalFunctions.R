load_expected <- function(mode) {
  mono <- readRDS(file.path("test-data", paste0(mode, "_recal_mono.rds")))
  iso <- readRDS(file.path("test-data", paste0(mode, "_recal_iso.rds")))
  recallist <- readRDS(file.path("test-data", paste0(mode, "_recal_recallist.rds")))
  return(list(Mono = mono, Iso = iso, RecalList = recallist))
}

update_expected <- function(actual, mode) {
  saveRDS(actual$Mono, file.path("test-data", paste0(mode, "_recal_mono.rds")))
  saveRDS(actual$Iso, file.path("test-data", paste0(mode, "_recal_iso.rds")))
  saveRDS(actual$RecalList, file.path("test-data", paste0(mode, "_recal_recallist.rds")))
}

patrick::with_parameters_test_that("Recal works",
  {
    peaks <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))
    recallist <- readRDS(file.path("test-data", paste0(mode, "_recallist.rds"))) |> dplyr::arrange_at("Series Score")

    actual <- MFAssignR::Recal(
      unambig,
      peaks = peaks$Mono,
      isopeaks = peaks$Iso,
      mzRange = 30,
      mode = mode,
      SN = 0,
      series1 = recallist$Series[1],
      series2 = recallist$Series[2],
      series3 = recallist$Series[3],
      series4 = recallist$Series[4],
      series5 = recallist$Series[5]
    )

    expected <- load_expected(mode)

    expect_equal(actual$Mono, expected$Mono)
    expect_equal(actual$Iso, expected$Iso)
    expect_equal(actual$RecalList, expected$RecalList)
  },
  mode = c("neg", "pos"),
)

patrick::with_parameters_test_that("RecalList works",
  {
    unambig <- readRDS(file.path("test-data", paste0(mode, "_cho_unambig.rds")))
    expected <- readRDS(file.path("test-data", paste0(mode, "_recallist.rds")))

    actual <- MFAssignR::RecalList(unambig)

    actual <- actual %>% dplyr::select(-"Series Index")
    expected <- expected %>% dplyr::select(-"Series Index")
    actual_sorted <- dplyr::arrange_at(actual, "Series")
    expected_sorted <- dplyr::arrange_at(expected, "Series")
    
    expect_equal(actual_sorted, expected_sorted)
  },
  mode = c("neg", "pos"),
)
