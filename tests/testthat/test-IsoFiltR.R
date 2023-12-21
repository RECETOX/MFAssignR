patrick::with_parameters_test_that("IsoFiltR works",
  {
    raw <- read.csv(file.path("test-data", paste0("QC1_1_", toupper(mode), "_500.csv")))
    expected <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    actual <- IsoFiltR(raw)
    expect_equal(actual, expected)
  },
  mode = c("pos", "neg")
)
