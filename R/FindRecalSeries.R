#' Filters the input dataframe
#' This function filters the input dataframe based on abundance score threshold and peak distance threshold; and 
#' computes the length of the series.
#' @param df Input dataframe = an output from RecalList, containing recalibrant CH2 series.
#' @param abundance_score_threshold A threshold for filtering abundance score parameter. The series with higher values #' are better. Default value is 100.
#' @param peak_distance_threshold A threshold for the peak distance parameter. The closer this value is to 1, the
#' better. 
#' @return A filtered dataframe.

filter_recal_series <- function(df, abundance_score_threshold, peak_distance_threshold) {
  df <- df %>%
    filter(Abundance.Score > abundance_score_threshold) %>%
    filter(Peak.Distance < peak_distance_threshold) %>%
    separate(col = Mass.Range, into = c('Min.Mass.Range', 'Max.Mass.Range'), sep = "-") %>%
    mutate(Min.Mass.Range = as.numeric(Min.Mass.Range), 
         Max.Mass.Range = as.numeric(Max.Mass.Range)) %>%
    mutate(Series.Length = Max.Mass.Range - Min.Mass.Range)

  return(df)
}

# Compute the scores
compute_scores <- function(combination) {
  series <- paste0(combination$Series)
  total_abundance <- sum(combination$Abundance.Score)
  total_series_length <- sum(combination$Series.Length)
  peak_score <- sum(1/(combination$Peak.Score))  
  peak_distance_proximity <- sum(1/(combination$Peak.Distance - 1))  
  coverage <- sum(combination$Max.Mass.Range - pmax(combination$Min.Mass.Range, lag(combination$Max.Mass.Range, default = 0)))
  coverage_percent <- coverage/((global_max+100) - (global_min-100))*100
  
  return(list(
    total_abundance = total_abundance,
    total_series_length = total_series_length,
    peak_proximity = peak_score,
    peak_distance_proximity = peak_distance_proximity,
    series = series,
    coverage = coverage,
    coverage_percent = coverage_percent
  ))
}

compute_coverage <- function(subset, global_max, global_min) {
  subset <- subset[order(subset$Min.Mass.Range), ]

  # Initialize the coverage and the end of the last segment
  total_coverage <- 0
  last_end <- -Inf

  # Iterate through the intervals
  for (i in 1:nrow(subset)) {
    current_start <- subset$Min.Mass.Range[i]
    current_end <- subset$Max.Mass.Range[i]
  
    if (current_start > last_end) {
      # Non-overlapping segment, add the full length to coverage
      total_coverage <- total_coverage + (current_end - current_start)
    } else {
      # Overlapping segment, add only the non-overlapping portion
      total_coverage <- total_coverage + max(0, current_end - last_end)
    }
  
    # Update the last_end to the current segment's end
    last_end <- max(last_end, current_end)
  }
  coverage_percent <- total_coverage/(global_max - global_min)*100
  return(coverage_percent)
}

covers_range <- function(subset, global_min, global_max, coverage_threshold) {
  coverage <- sum(subset$Max.Mass.Range - pmax(subset$Min.Mass.Range, lag(subset$Max.Mass.Range, default = 0)))
  coverage_percent <- coverage/((global_max-100) - (global_min+100))*100
  return(coverage_percent >= coverage_threshold)
}

find_series_combinations <- function(df, iter, global_min, global_max) {
  # Create empty list for scored combinations
  scores <- list()
  
  for (i in 1:nrow(iter)) {
    comb <- iter[i, ]
    subset <- df[comb, ]
    if (covers_range(subset, global_min, global_max, coverage_threshold)) {
        comb_score <- compute_scores(subset)
        scores <- append(scores, list(comb_score))
    }
   
}
}

#' Attempts to find most suitable series for recalibration.
#'
#' This function takes on input the CH2 homologous recalibration series, which are provided by the RecalList function #' and tries to find the most suitable series combination for recalibration based on the following criteria:
#' 1) Series should cover the full mass spectral range
#' 2) Series should be optimally long and combined have a “Tall Peak” at least every 100 m/z.
#' 3) Abundance score: the higher, the better
#' 4) Peak score: the closer to 0, the better
#' 5) Peak Distance: the closer to 1, the better
#' 6) Series Score: the closer to this value, the better 
#' 
#' The recal function can take up to 10 series - due to the size of the search space when looking for combinations of 10
#' elements, a pre-filtering is done: only the series which have Abundance score > 100 are considered and the one #'
#' having Peak Distance < 2. 
#' Combinations of 5 series are assembled, scores are computed for other metrics (in case of Peak proximity and Peak
#' distance, an inverted score is computed) and these are summed. Finally, top 10 unique series having the highest
#' score are outputted.
#' 
#' Warning: this step is in general computationally demanding, for ~30 series it took around 30 min. 
#' 
#' @param df An output from RecalList, containing recalibrant CH2 series.
#' @param abundance_score_threshold A threshold for filtering abundance score parameter. The series with higher values #' are better. Default value is 100.
#' @param peak_distance_threshold A threshold for the peak distance parameter. The closer this value is to 1, the
#' better.
#' @param tolerance A tolerance value to compute the global minimum and maximum. We expect that there is a low
#' probability that the true minimal/maximal m/z value of the dataset will be in a top-scoring series. Therefore we
#' need to set a reasonable tolerance, which will allow us to cover the most of the m/z range. Global minimum is then
#' computed as true minimum + tolerance ; global maximum as true maximum - tolerance. In case tolerance is set to 0,
#' true minimum and maximum will be used. Default value is 100.
#' @param combination_subset Combinations of how many series should be computed. Default is 5, Recal function can take
#' up to 10 series, but the more combinations, the longer computing time is expected (growing exponentially)
#' @param coverage_threshold How many % of the m/z range should be covered. Default is 90 %.
#' 
#' @return A dataframe of 10 best-scoring series.

find_series <- function(df, 
                       tolerance, 
                       combination_subset,
                       abundance_score_threshold,
                       peak_distance_threshold,
                       coverage_threshold) {

  # Arrange the data
  df <- filter_recal_series(df, abundance_score_threshold, peak_distance_threshold)

  # Compute the global minimum and maximum (range of a dataset)
  # We need to add some tolerance, because there is low chance full 100% would be covered
  global_min <- min(df$Min.Mass.Range) + tolerance
  global_max <- max(df$Max.Mass.Range) - tolerance

  # Create all combinations of ions
  iter <- combinations(nrow(df), combination_subset, v = 1:nrow(df))

  # Find the most appropriate combination of series
  scores_df <- find_series_combinations(df, iter, global_min, global_max, coverage_threshold)

# Iterate over combinations and score them
for (i in 1:nrow(coversRangeTrue)) {
  comb <- iter[i, ]
  subset <- df[comb, ]

}

# Append all scored combinations into a dataframe
scores_df <- do.call(rbind, lapply(scores, as.data.frame))

# Filter for the 10 top scoring series
finalSeries <- scores_df %>%
  filter(coverage_percent > 90) %>%
  rowwise() %>%
  mutate(sum_score = sum(total_abundance, total_series_length, peak_score, peak_distance_proximity, coverage_percent)) %>%
  arrange(desc(sum_score)) %>%
  filter(!duplicated(series)) %>%
  head(10)

# Return the top scoring series
return(finalSeries)
}