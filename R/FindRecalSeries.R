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
#' @return A dataframe of 10 best-scoring series.

findSeries <- function(df) {

# Arrange the data
df <- df %>%
  separate(col = Mass.Range, into = c('Min.Mass.Range', 'Max.Mass.Range'), sep = "-") %>%
  mutate(Min.Mass.Range = as.numeric(Min.Mass.Range), 
         Max.Mass.Range = as.numeric(Max.Mass.Range)) %>%
  mutate(Series.Length = Max.Mass.Range - Min.Mass.Range) %>%
  filter(Abundance.Score > 100) %>%
  filter(Peak.Distance < 2) 

# Compute the global minimum and maximum (range of a dataset)
# We need to add some tolerance, because there is low chance full 100% would be covered

tolerance <- 100
global_min <- min(df$Min.Mass.Range) + tolerance
global_max <- max(df$Max.Mass.Range) - tolerance

# Create all combinations of ions
iter <- combinations(nrow(df), 5, v = 1:nrow(df))

# Helper dataframe with information which combinations do cover range
coversRange <- data.frame(iter, coversRange = 0)

# Check if the combinations cover the whole data range
for (i in 1:nrow(iter)) {
  comb <- iter[i, ]
  subset <- df[comb, ]
  local_min <- min(subset$Min.Mass.Range)
  local_max <- max(subset$Max.Mass.Range)
  if (local_min <= global_min & local_max >= global_max) {
    coversRange$coversRange[i] <- 1
  } 
}

# Subset only those, which cover whole range
coversRangeTrue <- coversRange[coversRange$coversRange == 1, ]

# Compute the scores
score_combination <- function(combination) {
  series <- paste0(combination$Series)
  total_abundance <- sum(combination$Abundance.Score)
  total_series_length <- sum(combination$Series.Length)
  peak_proximity <- sum(1/(combination$Peak.Score))  
  peak_distance_proximity <- sum(1/(combination$Peak.Distance - 1))  
  coverage <- sum(combination$Max.Mass.Range - pmax(combination$Min.Mass.Range, lag(combination$Max.Mass.Range, default = 0)))
  coverage_percent <- coverage/((global_max+100) - (global_min-100))*100
  
  return(list(
    total_abundance = total_abundance,
    total_series_length = total_series_length,
    peak_proximity = peak_proximity,
    peak_distance_proximity = peak_distance_proximity,
    series = series,
    coverage = coverage,
    coverage_percent = coverage_percent
  ))
}

# Create empty list for scored combinations
scores <- list()

# Iterate over combinations and score them
for (i in 1:nrow(coversRangeTrue)) {
  comb <- iter[i, ]
  subset <- df[comb, ]
  comb_score <- score_combination(subset)
  scores <- append(scores, list(comb_score))
}

# Append all scored combinations into a dataframe
scores_df <- do.call(rbind, lapply(scores, as.data.frame))

# Filter for the 10 top scoring series
finalSeries <- scores_df %>%
  filter(coverage_percent > 90) %>%
  rowwise() %>%
  mutate(sum_score = sum(total_abundance, total_series_length, peak_proximity, peak_distance_proximity, coverage_percent)) %>%
  arrange(desc(sum_score)) %>%
  filter(!duplicated(series)) %>%
  head(10)

# Return the top scoring series
return(finalSeries)
}