patrick::with_parameters_test_that("IsoFiltR works",
  {
    raw <- read.csv(file.path("test-data", paste0("QC1_1_", toupper(mode), "_500.csv")))
    actual <- IsoFiltR(raw)
    expected <- readRDS(file.path("test-data", paste0(mode, "_iso.rds")))
    expect_equal(actual, expected)
  },
  mode = c("pos", "neg")
)

test_that("IsoFiltR works on recetox-aplcms output", {
  raw <- read.table(file.path("test-data", "21_qc_no_dil_milliq.txt"), header=TRUE, sep="\t")
  expected_path <- file.path("test-data", "21_qc_no_dil_milliq_iso.rds")

  actual <- IsoFiltR(raw)

  expected <- readRDS(expected_path)
  expect_equal(actual, expected)
})

patrick::with_parameters_test_that("filtered_data_is_empty works", {
  expect_equal(filtered_data_is_empty(test_data), expected)

  },
  patrick::cases(
    empty = list(test_data=list(content=c()), expected=TRUE),
    not_empty = list(test_data=list(content=c(1,2,3)), expected=FALSE)
  )
)