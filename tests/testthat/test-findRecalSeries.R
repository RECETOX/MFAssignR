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

test_that("Compute subsets from combinations", {
  df <- head(readRDS(file.path("test-data", "pos_recalSeries.rds")), 5)
  expected <- readRDS("test-data/combination_subsets.rds")
  actual <- compute_subsets(df, 3)
  expect_equal(actual, expected)
})