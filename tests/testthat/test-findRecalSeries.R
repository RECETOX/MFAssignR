test_that("Compute coverage works", {
  pos_recalList <- readRDS(file.path("test-data", "pos_recalSeries.rds"))
  actual <- compute_coverage(pos_recalList, 499.326, 111.044)
  expected <- 23.21534
  expect_equal(actual, expected, tolerance = 0.0005)
})

test_that("Compute combinations works", {
  df <- data.frame(series = c(1:5))
  expected <- readRDS("test-data/combinations_simple.rds")
  actual <- compute_combinations(df, 3)
  expect_equal(actual, expected)
})

test_that("Compute subsets from combinations works", {
  df <- head(readRDS(file.path("test-data", "pos_recalSeries.rds")), 5)
  expected <- readRDS("test-data/combination_subsets.rds")
  actual <- compute_subsets(df, 3)
  expect_equal(actual, expected)
})

test_that("Filtering of the subsets work", {
  df <- readRDS("test-data/combination_subsets.rds")
  expected <- readRDS("test-data/filtered_subsets.rds")
  actual <- filter_subsets_based_on_coverage(df, 80, 206, 117)
  expect_equal(actual, expected)
}
)

patrick::with_parameters_test_that("Selection of the final series works",
  {
    df <- readRDS("test-data/scores_df.rds")
    expected <- readRDS(file.path("test-data", paste0("final_series", mode, ".rds")))
    actual <- find_final_series(df, 3, mode)
    expect_equal(actual, expected)
  },
  mode = c(TRUE, FALSE)
)
