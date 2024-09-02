test_that("Compute coverage works", {
  pos_recalList <- readRDS(file.path("test-data", "pos_recalSeries.rds"))
  actual <- compute_coverage(pos_recalList, 499.326, 111.044)
  expected <- 23.5797
  expect_equal(actual, expected, tolerance = 0.5)
})