test_that("Filtering input dataframe works", {
  pos_recalList <- head(readRDS(file.path("test-data", "pos_recallist.rds")), 10)
  expected <- readRDS("test-data/filtered_recallist.rds")
  actual <- filter_recal_series(pos_recalList, abundance_score_threshold = 0, peak_distance_threshold = 2)
  
  expect_equal(actual, expected)
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

test_that("Compute coverage works", {
  pos_recalList <- readRDS(file.path("test-data", "pos_recalSeries.rds"))
  actual <- compute_coverage(pos_recalList, 499.326, 111.044)
  expected <- 23.21534
  expect_equal(actual, expected, tolerance = 0.0005)
})

test_that("Filtering of the subsets based on coverage works", {
  df <- readRDS("test-data/combination_subsets.rds")
  expected <- readRDS("test-data/filtered_subsets.rds")
  actual <- filter_subsets_based_on_coverage(df, 80, 206, 117)
  expect_equal(actual, expected)
})

test_that("Computing final scores work", {
  df <- head(readRDS("test-data/scores_df.rds"), 6)
  expected <- readRDS("test-data/final_scores.rds")
  actual <- compute_final_score(df)
  expect_equal(actual, expected)
})

test_that("Computing scores works", {
  df <- readRDS("test-data/filtered_subsets.rds")[[1]]
  expected <- readRDS("test-data/computed_scores.rds")
  actual <- compute_scores(df)
  expect_equal(actual, expected)
})

patrick::with_parameters_test_that("Selection of the final series works", {
  df <- readRDS("test-data/scores_df_full.rds")
  expected <- readRDS(file.path("test-data", paste0("final_series", mode, ".rds")))
  n <- 3

  actual <- find_final_series(df, n, mode)
  if (mode == TRUE) {
    expect_equal(nrow(actual), 10)
  } else {
    expect_equal(nrow(actual), n)
  }
  expect_equal(actual, expected)
},
  mode = c(TRUE, FALSE)
)

patrick::with_parameters_test_that("FindRecalSeries function works", {
  df <- readRDS("test-data/pos_recallist.rds")
  expected <- readRDS(file.path("test-data", paste0("findRecalSeries", mode, ".rds")))
  n <- 3

  actual <- FindRecalSeries(
    df,
    global_min = 100,
    global_max = 500,
    number_of_combinations = 3,
    abundance_score_threshold = 100,
    peak_distance_threshold = 2,
    coverage_threshold = 60,
    fill_series = mode)
  expect_equal(actual, expected)
},
  mode = c(TRUE, FALSE)
)

test_that("FindRecalSeriesSimple works", {
  df <- readRDS("test-data/pos_recallist.rds")
  actual <- FindRecalSeriesSimple(df)

  expected_path <- file.path("test-data", "expected_FindRecalSeriesSimple.rds")
  expected <- readRDS(expected_path)

  expect_equal(actual, expected)
})