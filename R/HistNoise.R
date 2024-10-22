#' Estimates the noise level
#'
#' For raw mass spectral data using the method
#' developed by Zhurov et al. Anal. Chem. (2014).
#' It uses the natural log of the intensity to create a histogram which then
#' is used to determine the noise level for the data being analyzed.
#' The data being put into this function should be the raw mass list
#' output from a mass spectrum, if the noise peaks have been removed
#' already, it will not work.
#'
#' The value reported in the console will be the estimated noise level
#' and the plot will demonstrate where the cut is being applied to the
#' data. The value in the console can then be multiplied by whatever
#' value is desired in order to reach the value to be used to cut the data.
#'
#' @param df - dataframe of intensity and ion mass, column 1 should be mass,
#' column 2 should be intensity
#' @param SN - numeric value for situations where a predefined noise value
#' is desired, default is 0
#' @param bin - numeric value determining the binwidth of the histogram,
#' default is 0.01
#'
#' @return List containing numeric noise level ('Noise')
#' and natural log intensity histogram ('Hist')
#'
#'   Histogram shows where the cut is being applied
#'
#' @export
HistNoise <- function(df, SN = 0, bin = 0.01) {
    df <- df[c(2, 1)]

    names(df)[2] <- "mass"
    names(df)[1] <- "intensity"

    Histo <- ggplot2::ggplot(df, ggplot2::aes(x = log(intensity))) + ggplot2::geom_histogram(binwidth = bin)
    Histo
    Data <- ggplot2::ggplot_build(Histo)$data  # Extracts the data from the histogram
    Count <- Data[[1]]$count  # Extracts the count data from the histogram
    LogInt <- Data[[1]]$x  # Extracts the x axis data from the histogram
    Freqdf <- data.frame(Count, LogInt)
    mode <- Freqdf[Count == max(Count), ]
    mode <- (mode[, 2])
    Valley <- Freqdf[(LogInt > mode & LogInt < (mode + 2)), ]  # Finds the valley between the noise max and the analyte
    ValMin <- Valley[Valley$Count == min(Valley$Count), ]  # Finds the minimum point of the valley
    ValMin$Ave <- mean(ValMin$LogInt)  # Average intensity of the minimum.
    Int <- ValMin[, 3]
    Int <- mean(Int)
    OutInt <- exp(Int)

    if (SN == 0) {
        df$Index <- "Bad"
        df$Index <- replace(df$Index, df$intensity > OutInt, "Good")
        print(OutInt)
        SNout <- OutInt
    }

    if (SN != 0) {
        df$Index <- "Bad"
        df$Index <- replace(df$Index, df$intensity > exp(SN), "Good")
        print(exp(SN))
        SNout <- exp(SN)
    }

    Freq2 <- ggplot2::ggplot(df, ggplot2::aes(x = log(intensity))) + ggplot2::geom_histogram(ggplot2::aes(y = ..count.., color = Index, fill = Index), binwidth = bin) + ggplot2::labs(x = "ln(Intensity)")

    Output <- list(Noise = SNout, Hist = Freq2)
    Output
}
