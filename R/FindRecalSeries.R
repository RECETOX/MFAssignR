#' Filters the input dataframe
#' This function filters the input dataframe based on Abundance.Score threshold and peak distance threshold; and
#' computes the length of the series.
#' @param df DataFrame An output from RecalList, containing recalibrant CH2 series.
#' @param abundance_score_threshold Float A threshold for filtering Abundance.Score parameter. The series with higher values #' are better. Default value is 100.
#' @param peak_distance_threshold Float A threshold for the peak distance parameter. The closer this value is to 1, the
#' better.
#' @return DataFrame A filtered dataframe.
#' @importFrom dplyr %>%
#' @export
filter_recal_series <- function(df, abundance_score_threshold, peak_distance_threshold) {
  df <- df %>%
    dplyr::filter(Abundance.Score > abundance_score_threshold) %>%
    dplyr::filter(Peak.Distance < peak_distance_threshold) %>%
    tidyr::separate(col = Mass.Range, into = c("Min.Mass.Range", "Max.Mass.Range"), sep = "-") %>%
    dplyr::mutate(
      Min.Mass.Range = as.numeric(Min.Mass.Range),
      Max.Mass.Range = as.numeric(Max.Mass.Range)
    ) %>%
    dplyr::mutate(Series.Length = Max.Mass.Range - Min.Mass.Range)

  return(df)
}

#' Computes all possible combinations of series
#' This function computes all possible combinations of series. The number of combinations is specified by user.
#' @param df DataFrame A pre-filtered dataframe containing series.
#' @param n Integer Number of combinations to compute.
#'
#' @return DataFrame A dataframe of row-index based combinations.
#' @export
compute_combinations <- function(df, n) {
  combs <- gtools::combinations(nrow(df), n, v = seq_len(nrow(df)))
  return(combs)
}

#' Computes subsets based on combinations.
#' This function selects subsets of the input data based on the combinations of series.
#' @param df DataFrame A pre-filtered dataframe containing series.
#' @param n Integer Number of combinations to compute.
#'
#' @return A list of subsets based on series combinations.
#' @export
compute_subsets <- function(df, n) {
  subsets <- apply(compute_combinations(df, n), 1, function(x) {
    df[x, ]
  })
  return(subsets)
}

#' Computes coverage of a subset.
#' This function computes m/z range coverage of a particular subset. The higher coverage, the better.
#' @param subset DataFrame A subset of pre-filtered input data, has a size of combination number.
#' @param global_min Float A lower bound of the instrument m/z range.
#' @param global_max Float A higher bound of the instrument m/z range.
#'
#' @return Coverage (in %) of a particular subset.
#' @importFrom dplyr %>%
#' @export
compute_coverage <- function(subset, global_max, global_min) {
  subset <- subset %>%
    dplyr::arrange(Min.Mass.Range)
  
  #subset <- subset[order(subset$Min.Mass.Range), ]

  # Initialize the coverage and the end of the last segment
  total_coverage <- 0
  last_end <- -Inf

  # Iterate through the intervals
  for (i in seq_len(nrow(subset))) {
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
  coverage_percent <- total_coverage / (global_max - global_min) * 100
  return(coverage_percent)
}

#' Filters subsets based on the coverage.
#' This function filters the subsets based on the coverage threshold.
#' @param subsets DataFrame A subset of pre-filtered input data, has a size of combination number.
#' @param coverage_threshold Integer How many % of the m/z range should be covered. Default is 90 %.
#' @param global_min Float A lower bound of the instrument m/z range.
#' @param global_max Float A higher bound of the instrument m/z range.
#'
#' @return List A list of subsets which pass the coverage threshold (their coverage is higher than the threshold)
#' @export
filter_subsets_based_on_coverage <- function(subsets, coverage_threshold, global_max, global_min) {
  filtered_subsets <- Filter(function(x) compute_coverage(x, global_max, global_min) > coverage_threshold, subsets)
  return(filtered_subsets)
}

#' Computes the scores for individual subsets.
#' This function computes individual scores for the particular combination of series.
#' @param combination DataFrame A subset of pre-filtered input data, has a size of combination number.
#'
#' @return List List of scores for individual parameters.
#' @export
compute_scores <- function(combination) {
  series <- paste0(combination$Series)
  total_abundance <- sum(combination$Abundance.Score)
  total_series_length <- sum(combination$Series.Length)
  peak_score <- sum(1 / (combination$Peak.Score))
  peak_distance_proximity <- sum(1 / (combination$Peak.Distance - 1))
  series_id <- paste(combination$Series, collapse = " ")

  return(list(
    series = series,
    total_abundance = total_abundance,
    total_series_length = total_series_length,
    peak_proximity = peak_score,
    peak_distance_proximity = peak_distance_proximity,
    series_id = series_id
  ))
}

#' Compute a sum score for the whole combination.
#' This function computes a score summing all parameter scores and arranges them in descending manner.
#' @param scores_df DataFrame A dataframe of scores.
#'
#' @return DataFrame A sorted dataframe containing the final score.
#' @importFrom dplyr %>%
#' @export
compute_final_score <- function(scores_df) {
  final_score <- scores_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sum_score = sum(total_abundance, total_series_length, peak_proximity, peak_distance_proximity)) %>%
    dplyr::arrange(desc(sum_score)) %>%
    dplyr::ungroup() 

  return(final_score)
}

#' Selects final series.
#' This function selects final best-scoring series. An user can choose, whether only a number of series corresponding
#' to the size of the combination will be returned, or 10 best-scoring 10 series will be returned.
#' @param scores_df DataFrame A dataframe of scores.
#' @param number_of_combinations Integer Number of combinations to compute.
#' @param fill_series Boolean If TRUE, top 10 unique series will be returned, otherwise only the series from the best
#' combination will be returned.
#'
#' @return DataFrame A dataframe of best-scoring series.
#' @importFrom dplyr %>%
#' @export
find_final_series <- function(scores_df, number_of_combinations, fill_series) {
  final_series <- compute_final_score(scores_df)

  if (fill_series == FALSE) {
    final_series <- final_series %>%
      dplyr::slice_head(n = number_of_combinations)
  } else {
    final_series <- final_series %>%
      dplyr::distinct(series, .keep_all = TRUE) %>%
      dplyr::slice_head(n = 10)
  }

  return(final_series)
}

#' Attempts to find most suitable series for recalibration.
#'
#' This function takes on input the CH2 homologous recalibration series, which are provided by the RecalList function #' and tries to find the most suitable series combination for recalibration based on the following criteria:
#' 1) Series should cover the full mass spectral range,
#' 2) Series should be optimally long and combined have a “Tall Peak” at least every 100 m/z,
#' 3) Abundance.Score: the higher, the better,
#' 4) Peak score: the closer to 0, the better,
#' 5) Peak Distance: the closer to 1, the better,
#' 6) Series Score: the closer to this value, the better.
#'
#' The recal function can take up to 10 series - due to the size of the search space when looking for combinations of 10
#' elements, a pre-filtering is done: only the series which have Abundance.Score > 100 are considered and the one #'
#' having Peak Distance < 2.
#' Combinations of 5 series are assembled, scores are computed for other metrics (in case of Peak proximity and Peak
#' distance, an inverted score is computed) and these are summed. Finally, top 10 unique series having the highest
#' score are outputted.
#'
#' Warning: this step is in general computationally demanding, for ~30 series it took around 30 min.
#'
#' @param df DataFrame An output from RecalList, containing recalibrant CH2 series.
#' @param global_min Float A lower bound of the instrument m/z range.
#' @param global_max Float A higher bound of the instrument m/z range.
#' @param number_of_combinations Integer Combinations of how many series should be computed. Default is 5, Recal function can
#' take up to 10 series, but the more combinations, the longer computing time is expected (growing exponentially)
#' @param abundance_score_threshold Float A threshold for filtering Abundance.Score parameter. The series with higher values #' are better. Default value is 100.
#' @param peak_distance_threshold Float A threshold for the peak distance parameter. The closer this value is to 1, the
#' better.
#' @param coverage_threshold Integer How many % of the m/z range should be covered. Default is 90 %.
#' @param fill_series Boolean If TRUE, top 10 unique series will be returned, otherwise only the series from the best
#' combination will be returned.
#'
#' @return A dataframe of n-10 best-scoring series.
#' @export
FindRecalSeries <- function(df,
                        global_min,
                        global_max,
                        number_of_combinations,
                        abundance_score_threshold,
                        peak_distance_threshold,
                        coverage_threshold,
                        fill_series) {
  # Pre-filter and arrange the data
  df <- filter_recal_series(df, abundance_score_threshold, peak_distance_threshold)

  # Create all combinations of ions
  subsets <- compute_subsets(df, number_of_combinations)

  # Filter the subsets based on coverage threshold
  subsets_filtered <- filter_subsets_based_on_coverage(subsets, coverage_threshold, global_max, global_min)

  # Compute the scores
  scores <- lapply(subsets_filtered, function(x) {
    compute_scores(x)
  })

  # Append all scored combinations into a dataframe
  scores_df <- do.call(rbind, lapply(scores, as.data.frame))

  # Filter for final series
  final_series <- find_final_series(scores_df, number_of_combinations, fill_series)

  # Return the top scoring series
  return(final_series)
}

#' Simple rewritten version of the FindRecalSeries function.
#' 
#' @description This function is not based on combinations of series but simply computes the scores and returns
#' the 10 best seires.
#' @param Recal data.frame A dataframe containing the various recal series.
#' @return A dataframe of n-10 best-scoring series.
#' @export
FindRecalSeriesSimple <- function(Recal) {
  Cal_Pick <- dplyr::filter(Recal, `Series.Score`>= 1 & `Peak.Distance` <= 3.3 & `Peak.Distance` >= 1)
  Cal_Pick <- dplyr::mutate(Cal_Pick, `Peak.Distance` = floor(`Peak.Distance`), `Series.Score` = round(`Series.Score`, 1))

  #The weighting for each of the terms was determined experimentally and via my experience choosing series and what is most important.
  Cal_Pick <- dplyr::mutate(
    Cal_Pick,
    Num_Weight = `Number.Observed` / max(`Number.Observed`) *25,
    Abund_weight = (((`Abundance.Score`) / max(`Abundance.Score`))*20),
    TallP_Weight = `Tall.Peak`/30,
    PeakS_Weight = abs((`Peak.Score`* 5 - 10)),
    PeakD_Weight = -3 * `Peak.Distance` + 10,
    Series_Weight = (`Series.Score` / (( `Series.Score` - 0.5 )))*5
  )

  Cal_Pick <- dplyr::mutate(Cal_Pick, Total_score = Num_Weight + PeakS_Weight + PeakD_Weight + Abund_weight+ Series_Weight + TallP_Weight)  
  Cal_Pick <- Cal_Pick[order(-Cal_Pick$Total_score),]

  Cal_Pick_HM <- tidyr::separate(Cal_Pick, `Mass.Range`, into = c("Low", "High"), sep = "-", remove = FALSE)  
  Top10 <- dplyr::slice(Cal_Pick_HM, c(1:10))

  Cal_Pick_HM2 <- dplyr::filter(Cal_Pick_HM, Low > as.numeric(max(Top10$High))-50 & Low < as.numeric(max(Top10$High))-10)  
  HM_Cal <- dplyr::filter(Cal_Pick_HM2, Total_score == max(Total_score))

  Cal_Pick_LM <- dplyr::filter(Cal_Pick_HM, (Low > as.numeric(min(Cal_Pick_HM$Low)) & Low < as.numeric(min(Cal_Pick_HM$Low))+20) & High > as.numeric(min(Top10$Low)))
  LM_Cal <- dplyr::filter(Cal_Pick_LM, Total_score == max(Total_score))

  Picked_Cal_series <- dplyr::bind_rows(LM_Cal, Top10, HM_Cal)
  Picked_Cal_series <- dplyr::select(Picked_Cal_series, -c(5,6))

  Picked_Cal_series
}