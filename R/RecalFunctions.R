#' Plot Spectrum Function
#'
#' This function generates a mass spectrum plot using ggplot2.
#'
#' @param plotpeak A data frame containing information for the first set of segments.
#' @param RecalPlot A data frame containing information for the second set of segments.
#'
#' @return A ggplot object representing the mass spectrum plot.
plot_spectrum <- function(plotpeak, RecalPlot) {
  ggplot() +
    geom_segment(data = plotpeak, size = 0.7, aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance"), alpha = 0.3, color = "grey57") +
    geom_segment(data = RecalPlot, size = 1.2, aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "Abundance"), color = "blue") +
    theme_bw() +
    labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "Series") +
    theme(
      axis.title = element_text(size = 15, face = "bold"),
      strip.text = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(face = "bold", size = 15),
      legend.text = element_text(face = "bold", size = 15),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_blank(),
      plot.title = element_text(size = 16, face = "bold")
    )
}

#' Process Mass List
#'
#' This function processes a data frame containing mass-related information.
#'
#' @param data A data frame containing at least the "Exp_mass" column.
#' @return The input data frame with additional calculated columns.
#'
#' @details
#' The function calculates the following columns:
#'   - NM: Rounded Exp_mass
#'   - KM_O: Exp_mass * (16 / 15.9949146223)
#'   - KMD_O: Rounded NM - KM_O
#'   - z_O: Rounded Exp_mass %% 16 - 16
#'   - KM_H2: Exp_mass * (2 / 2.01565)
#'   - KMD_H2: Rounded NM - KM_H2
#'   - z_H2: Rounded Exp_mass %% 2 - 2
#'
#' @examples
#' data <- data.frame(Exp_mass = c(100, 150, 200))
#' processed_data <- processMassList(data)
#' @export
processMassList <- function(data) {
  data$NM <- round(data$Exp_mass)

  # Process O
  data$KM_O <- data$Exp_mass * (16 / 15.9949146223)
  data$KMD_O <- round(data$NM - data$KM_O, 3)
  data$z_O <- round(data$Exp_mass) %% 16 - 16

  # Process H2
  data$KM_H2 <- data$Exp_mass * (2 / 2.01565)
  data$KMD_H2 <- round(data$NM - data$KM_H2, 3)
  data$z_H2 <- round(data$Exp_mass) %% 2 - 2

  return(data)
}

#' Process Known Compounds
#'
#' This function processes known compounds by merging them with a rest dataset and
#' updating abundance and type columns based on specified columns and limits.
#'
#' @param rest The dataset containing the rest of the compounds.
#' @param known The dataset containing known compounds.
#' @param kmd_col The column name for KMD in both datasets.
#' @param z_col The column name for Z in both datasets.
#' @param num_col The column name to store the calculated number of compounds.
#' @param element_mass The mass of either O16 or H2.
#' @param type The column name to store the type of compound.
#' @param step_limit The limit for considering compounds.
#' @param remove_indices Indices of rows to be removed from the result.
#'
#' @return A processed dataset with updated abundance and type columns.
processKnown <- function(rest, known, kmd_col, z_col, num_col, element_mass, type, step_limit, remove_indices) {
  names(known)[2] <- "base_mass"

  step_result <- merge(rest, known, by.x = c(kmd_col, z_col), by.y = c(kmd_col, z_col))
  step_result[[num_col]] <- round((step_result$Exp_mass - step_result$base_mass) / element_mass)
  step_result[[type]] <- step_result[[type]] + step_result[[num_col]]
  step_result$Type <- type
  step_result$form <- paste(step_result$C, step_result$H, step_result$O, step_result$N, step_result$S,  step_result$P, step_result$E, step_result$S34, step_result$N15, step_result$D, step_result$Cl, step_result$Fl, step_result$Cl37, step_result$M, step_result$NH4, step_result$POE, step_result$NOE, sep = "_")
  step_result <- step_result[abs(step_result[[num_col]]) <= step_limit, ]
  step_result <- step_result[-remove_indices]
  return(step_result)
}

#' Generates a plot to visualize recalibrant series, recalibrates two mass lists,
#' and produces a list of the chosen recalibrants.
#'
#'
#' Recal() takes the data frame output from MFAssign functions, the outputs of
#' IsoFiltR(), and the chosen recalibrant series to generate a plot for the qualitative
#' assessment of the recalibrants, recalibrate the "Mono" and "Iso" outputs from IsoFiltR,
#' and prepare a data frame containing the chosen recalibrants.
#'
#'
#' This function can handle up to 10 homologous series, though it will work with fewer.
#' It is important for recalibrant masses to cover the entire mass range of interest,
#' and they should be among the most abundant peaks in their region of the spectrum,
#' this function helps with visualizing the recalibrants.
#'
#' This function also recalibrates up to two mass lists using the chosen recalibrants.
#' It is best to use the "Mono" and "Iso" outputs of the IsoFiltR() function.
#'
#' This recalibration method uses the first step from the paper from Kozhinov et al.
#' 2014 (Anal. Chem.) to estimate the correction term, and uses the idea of a segmented
#' or "walking" recalibration from Savory et al. 2011 (Anal. Chem) to further improve
#' the recalibration.
#'
#'
#' These parameters will be used to identify more potential recalibrants, and the
#' tallest peaks within a user defined mass range will be used as recalibrants.
#'
#' @param df data frame output from MFAssign() or MFAssignCHO()
#'
#' @param peaks data frame with two columns for recalibration, generally contains "Mono" output from IsoFiltR()
#'
#' @param isopeaks data frame with two columns, generally containing the "Iso" output from IsoFiltR()
#'
#' @param mode character "neg" or "pos" depending on the ionization mode of the data
#'
#' @param SN numeric value of the signal-to-noise cut point for making the plot, default is 0
#'
#' @param mzRange numeric value defining the mass windows used for the segmented recalibration, default is 30
#'
#' @param series(1-10) character the recalibrant series, "O7_Na_4" for example, default is NA
#'
#' @param step_O numeric value defining the number of oxygen "steps" for formula extension, default is 3
#'
#' @param step_H2 numeric value defining the number of H2 "steps" for formula extension, default is 5
#'
#' @param CalPeak numeric value defining the maximum allowed recalibrant peaks per mzRange defined segment,
#' default is 150
#'
#'
#' @return list(MZ, Mono, Iso, RecalOut) contains the mass spectrum (MZ), recalibrated "Mono" data frame,
#'  recalibrated "Iso" data frame, and the recalibrants list (RecalOut)
#'
#' @export
Recal <- function(df,
                  peaks,
                  isopeaks = "none",
                  mode,
                  SN = 0,
                  mzRange = 30,
                  series1 = NA, series2 = NA,
                  series3 = NA, series4 = NA,
                  series5 = NA, series6 = NA,
                  series7 = NA, series8 = NA,
                  series9 = NA, series10 = NA,
                  step_O = 3,
                  step_H2 = 5,
                  CalPeak = 150) {
  cols <- ncol(peaks)
  min <- round(min(peaks[1]))
  max <- round(max(peaks[1]))

  RecalList <- data.frame(
    series = c(
      series1, series2, series3, series4, series5, series6, series7, series8, series9,
      series10
    ),
    Tag = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10")
  )
  RecalList <- RecalList[!is.na(RecalList$series), ]

  df$Adduct <- "H"
  df$Adduct <- replace(df$Adduct, df$M > 0, "Na")
  df$Adduct <- replace(df$Adduct, df$POE == 1, "OE")
  df$Adduct <- replace(df$Adduct, df$NOE == 1, "OE")
  df$series <- paste(df$class, df$Adduct, df$DBE, sep = "_")
  df$mode <- mode
  names(df)[1:2] <- c("Abundance", "Exp_mass")

  if (cols == 2) {
    isopeaks <- if(typeof(isopeaks) == "character"){data.frame(mass = 1, Abundance = 1, Tag = "X")}else{isopeaks}
    # isopeaks <- if(isopeaks == "none") data.frame(mass = 1, Abundance = 1, Tag = "X") else isopeaks
    isopeaks <- data.frame(isopeaks)
    isopeaks <- isopeaks[, c(1:3)]
    isopeaks <- setNames(isopeaks, c("mass", "Abundance", "Tag"))
    peaks <- setNames(peaks, c("Abundance", "mass"))

    isopeaks$Tag2 <- "Iso"
    peaks$Tag <- "X"
    peaks$Tag2 <- "Mono"

    # peaks <- if(isopeaks != "none") rbind(peaks, isopeaks) else peaks
    peaks <- ifelse(isopeaks != "none", rbind(peaks, isopeaks), peaks) # Change 12/6/19
    peaks <- data.frame(peaks)
    peaks <- peaks[, c(1:4)]
    peaks <- setNames(peaks[c(2, 1, 3, 4)], c("Abundance", "mass", "Tag", "Tag2"))
    isopeaks <- isopeaks[c(2, 1, 3, 4)]

    # Merges recalibrant list to df in order to determine which recalibrants are in the data frame.
    RecalList2 <- merge(df, RecalList, by.x = "series", by.y = "series")


    # cehck <- Iso %>% filter(!is.na(exp_mass))
    ##########################################
    # Monoisotopic Peaks
    # Prepares the recalibrant masses for use in recalibration steps.
    RecalList <- RecalList2[c(
      "Abundance", "Exp_mass", "C", "H", "O", "N", "S", "P", "E",
      "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z", "C13_mass"
    )]
    RecalList <- processMassList(RecalList)

    Rest <- df[c("Abundance", "Exp_mass")]
    Rest <- processMassList(Rest)

    ############
    # Picking recalibrants with series
    # Process known_O
    O_mass <- 15.9949146223
    Step2 <- processKnown(Rest, RecalList[c(1:21, 24, 25)], "KMD_O", "z_O", "O_num", O_mass, "O", step_O, c(10, 31))

    # Process known_H2
    H_mass <- 2.01565
    Step3 <- processKnown(Rest, RecalList[c(1:21, 27, 28)], "KMD_H2", "z_H2", "H2_num", H_mass, "H2", step_H2, c(10, 31))

    Out <- rbind(Step2, Step3)
    Out2 <- dplyr::distinct(Out, Exp_mass)

    NewRecal <- merge(df, Out2, by = "Exp_mass")
    NewRecal <- NewRecal[c("Abundance", "Exp_mass", "formula", "theor_mass", "C13_mass", "C13_abund")]
    # Good to here
    #####################
    # Picking recalibrants by abundance
    NewRecal$Bin <- cut(NewRecal$Exp_mass, breaks = seq(min, max, by = mzRange))

    FinalRecal2 <- dplyr::group_by(NewRecal, Bin)
    # FinalRecal2 <- dplyr::top_n(FinalRecal, obs, Abundance)

    NewRecal_mono <- FinalRecal2[c(1:4)] # Changed 5/13/20
    NewRecal_iso <- FinalRecal2[c(5, 3, 4, 6)] # Changed 5/13/20
    NewRecal_iso <- NewRecal_iso[NewRecal_iso$C13_mass > 0, ]
    names(NewRecal_iso)[1] <- "Exp_mass"
    names(NewRecal_iso)[4] <- "Abundance" # Changed 5/13/20
    NewRecal_iso$theor_mass <- NewRecal_iso$theor_mass + 1.003355
    NewRecal <- rbind(NewRecal_mono, NewRecal_iso)
  }
  ######################################################
  # Version for LC-MS Data Lists
  if (cols == 3) {
    isopeaks <- if(typeof(isopeaks) == "character"){data.frame(mass = 1, Abundance = 1, RT = 0, Tag = "X")}else{isopeaks}

    isopeaks <- data.frame(isopeaks)
    isopeaks <- isopeaks[, c(1:4)]
    isopeaks <- setNames(isopeaks, c("mass", "Abundance", "RT", "Tag"))
    peaks <- setNames(peaks, c("mass", "Abundance", "RT"))

    isopeaks$Tag2 <- "Iso"
    peaks$Tag <- "X"
    peaks$Tag2 <- "Mono"


    peaks <- ifelse(isopeaks != "none", rbind(peaks, isopeaks), peaks) # Change 12/6/19
    peaks <- data.frame(peaks)
    peaks <- peaks[, c(1:5)]
    peaks <- setNames(peaks, c("mass", "Abundance", "RT", "Tag", "Tag2"))


    peaks <- peaks[c(2, 1, 3, 4, 5)]
    # isodummy <- data.frame
    isopeaks <- isopeaks[c(2, 1, 3, 4, 5)]

    # Merges recalibrant list to df in order to determine which recalibrants are in the data frame.
    RecalList2 <- merge(df, RecalList, by.x = "series", by.y = "series")

    ##########################################
    # Monoisotopic Peaks
    # Prepares the recalibrant masses for use in recalibration steps.
    RecalList <- RecalList2[c(
      "Abundance", "Exp_mass", "RT", "C", "H", "O", "N", "S", "P", "E",
      "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z", "C13_mass"
    )]
    RecalList <- processMassList(RecalList)

    Rest <- df[c("Abundance", "Exp_mass", "RT")]
    # Assuming RecalList is your data frame
    Rest <- processMassList(Rest)

    ############
    # Picking recalibrants with series

    knownO <- RecalList[c(1:22, 25, 26)]
    names(knownO)[2] <- "base_mass"
    Step2 <- merge(Rest, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
    Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass)) / 15.9949146223)
    Step2$O <- Step2$O + Step2$O_num
    Step2$Type <- "O"
    Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
      Step2$N15, Step2$D, Step2$Cl, Step2$Fl, Step2$Cl37, Step2$M, Step2$NH4,
      Step2$POE, Step2$NOE,
      sep = "_"
    )
    Step2 <- Step2[abs(Step2$O_num) <= step_O, ]
    Step2 <- Step2[-c(11, 13, 33)] #

    knownH2 <- RecalList[c(1:22, 28, 29)]
    names(knownH2)[2] <- "base_mass"
    Step3 <- merge(Rest, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
    Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass)) / 2.01565)
    Step3$H <- Step3$H + 2 * Step3$H2_num
    Step3$Type <- "H2"
    Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
      Step3$N15, Step3$D, Step3$Cl, Step3$Fl, Step3$Cl37, Step3$M, Step3$NH4,
      Step3$POE, Step3$NOE,
      sep = "_"
    )
    Step3 <- Step3[abs(Step3$H2_num) <= step_H2, ]
    Step3 <- Step3[-c(11, 13, 33)]

    Out <- rbind(Step2, Step3)
    Out2 <- dplyr::distinct(Out, Exp_mass, .keep_all = TRUE)
    Out2 <- Out2[c(3, 4)]
    names(Out2)[1] <- "Abundance"

    NewRecal <- merge(df, Out2, by = c("Exp_mass", "Abundance"))
    NewRecal <- NewRecal[c("Abundance", "Exp_mass", "RT", "formula", "theor_mass", "C13_mass", "C13_abund")]
    # Good to here
    #####################
    # Picking recalibrants by abundance
    NewRecal$Bin <- cut(NewRecal$Exp_mass, breaks = seq(min, max, by = mzRange))

    FinalRecal2 <- dplyr::group_by(NewRecal, Bin)


    NewRecal_mono <- FinalRecal2[c(1:5)] # Changed 5/13/20
    NewRecal_iso <- FinalRecal2[c(6, 3, 4, 5, 7)] # Changed 5/13/20
    NewRecal_iso <- NewRecal_iso[NewRecal_iso$C13_mass > 0, ]
    names(NewRecal_iso)[1] <- "Exp_mass"
    names(NewRecal_iso)[5] <- "Abundance" # Changed 5/13/20
    NewRecal_iso$theor_mass <- NewRecal_iso$theor_mass + 1.003355
    NewRecal <- rbind(NewRecal_mono, NewRecal_iso)
  }

  ########################################################
  # Setting up the recalibrant and mass lists so they can be recalibrated by section.
  # mzRange <- 50
  peaks <- peaks[!is.na(peaks$mass), ]
  Dummy <- data.frame(Values = seq(from = floor(min(peaks$mass)), to = ceiling(max(peaks$mass)), by = mzRange))
  Dummy2 <- nrow(Dummy)
  NewRecal$Range <- 1

  for (i in 1:Dummy2) { # Sets up the bins for the recalibrants
    I <- i

    NewRecal$Range <- replace(NewRecal$Range, NewRecal$Exp_mass >= min(NewRecal$Exp_mass) + (I - 1) * mzRange &
      NewRecal$Exp_mass < min(NewRecal$Exp_mass) + I * mzRange, I)
  }

  peaks$Range <- 1

  for (i in 1:Dummy2) { # Sets up the bins for the raw data
    I <- i

    peaks$Range <- replace(peaks$Range, peaks$mass >= min(peaks$mass) + (I - 1) * mzRange &
      peaks$mass < min(peaks$mass) + I * mzRange, I)
  }

  # This section calculates how many recalibrants are in each segment
  # test <- NewRecal
  NewRecal$Test <- 1
  for (i in 1:Dummy2) {
    I <- i

    NewRecal$Test <- replace(NewRecal$Test, NewRecal$Range == I, nrow(NewRecal[NewRecal$Range == as.numeric(I), ]))
  }


  ####################################### 33
  # Fixing the problem with too many recalibrants in a single section. 5/13/20
  ## THIS IS WHERE I AM NOW

  Orig_Range <- data.frame(Range = unique(peaks$Range)) # NEW
  # Gap <- NewRecal[NewRecal$Test < 3,]  #NEW
  NewRecal <- NewRecal[NewRecal$Test > 2, ]
  # Gap <- Gap[Gap$Range < max(NewRecal$Range),]
  Step_List <- unique(NewRecal[c(5)])
  if (cols == 3) {
    Step_List <- unique(NewRecal[c(6)])
  }

  Step_List$Tag <- "Good"

  GapCheck <- merge(Orig_Range, Step_List, by.x = "Range", by.y = "Range", all = T)
  GapCheck <- GapCheck[is.na(GapCheck$Tag) & GapCheck$Range < max(NewRecal$Range), ]

  Range_List <- split(NewRecal, f = NewRecal$Range)

  i <- 1
  Range_ListF <- data.frame(Abundance = 0, Exp_mass = -1, formula = "X", theor_mass = -1, Range = 0, Test = 1)
  if (cols == 3) {
    Range_ListF <- data.frame(Abundance = 0, Exp_mass = -1, RT = -42, formula = "X", theor_mass = -1, Range = 0, Test = 1)
  }

  # Range_ListF <- list()
  for (i in 1:nrow(Step_List)) {
    Range_List2 <- Range_List[[i]]
    Range_List2 <- dplyr::top_n(Range_List2, CalPeak, Range_List2$Abundance)
    # Range_ListF[[length(Range_ListF)+1]] = Range_List2
    Range_ListF <- rbind(Range_ListF, Range_List2)
  }
  NewRecal <- Range_ListF[Range_ListF$Exp_mass > 0, ]
  # NewRecal <- NewRecal[-1]

  #######################################
  # Error message if the problem is having too many recalibrants
  ## START FROM HERE

  if (nrow(GapCheck) > 0) stop("Gap in recalibrant coverage, try adding more recalibrant series") # New

  NewRecalx <- NewRecal[c(1:5)]
  if (cols == 3) {
    NewRecalx <- NewRecal[c(1:6)]
  } # LCMS

  NewRecal <- NewRecalx # LCMS

  # This is just setting up a dataframe for final export
  Align <- NewRecal
  names(Align)[2] <- "Exp_mass"
  FinalRecal2 <- merge(FinalRecal2, Align, by.x = c("Exp_mass", "Abundance"), by.y = c("Exp_mass", "Abundance")) # LCMS
  FinalRecal2x <- FinalRecal2[c(1:6)] # LCMS
  if (cols == 3) {
    FinalRecal2x <- FinalRecal2[c(1:7)]
  } # LCMS
  FinalRecal2 <- FinalRecal2x # LCMS

  if (cols == 2) {
    names(FinalRecal2)[3] <- "formula"
  } # LCMS
  if (cols == 3) {
    names(FinalRecal2)[4] <- "formula"
    names(FinalRecal2)[3] <- "RT"
  } # LCMS
  #########################
  # Calculating the weights needed for recalibration and calculating the recalibrated masses.
  if (cols == 2) {
    colnames(NewRecal)[c(2, 3, 4)] <- c("E_mass", "formula", "Th_mass")
    peaks_dum <- data.frame(Abundance = -2, mass = -2, Tag = "X", Tag2 = "Y", Range = -1)
    peaks_out <- peaks_dum
    Ej_dum <- data.frame(Ejsum = -2, Range = -2)
    Ej_out <- Ej_dum
    Dummy3 <- max(NewRecal$Range)
  }


  if (cols == 3) {
    colnames(NewRecal)[c(2, 4, 5)] <- c("E_mass", "formula", "Th_mass")
    peaks_dum <- data.frame(Abundance = -2, mass = -2, RT = -42, Tag = "X", Tag2 = "Y", Range = -1)
    peaks_out <- peaks_dum
    Ej_dum <- data.frame(Ejsum = -2, Range = -2)
    Ej_out <- Ej_dum
    Dummy3 <- as.numeric(nrow(data.frame(Value = 1:max(NewRecal$Range))))
  }
  ############################
  peakshigh <- peaks[peaks$Range > max(NewRecal$Range), ]
  peaks <- peaks[peaks$Range <= max(NewRecal$Range), ]

  # library(dplyr)
  check <- NewRecal[!duplicated(NewRecal$Range), ]
  check2 <- peaks[!duplicated(peaks$Range), ]

  if (nrow(check) != nrow(check2)) stop("Increase the value of mzRange")

  # j <- 1
  for (j in 1:(Dummy3)) {
    J <- j
    while (J <= Dummy3) {
      NewRecal_2 <- NewRecal[NewRecal$Range == J, ]
      peaks_2 <- peaks[peaks$Range == J, ]
      peaks_final <- rbind(peaks_dum, peaks_out)
      Ej_final <- rbind(Ej_dum, Ej_out)
      counter <- 1
      J <- J + Dummy3
      while (counter == 1) {
        counter <- counter + 1
        NewRecal_2$num <- 0:(nrow(NewRecal_2) - 1)
        NewRecal_2$weight <- factorial(nrow(NewRecal_2) - 1) /
          (factorial(NewRecal_2$num) * factorial((nrow(NewRecal_2) - 1) - NewRecal_2$num))

        NewRecal_2$mzweight <- NewRecal_2$weight * NewRecal_2$E_mass
        NewRecal_2$Ejweight <- NewRecal_2$weight * (NewRecal_2$E_mass - NewRecal_2$Th_mass) / NewRecal_2$Th_mass
        NewRecal_2$Ejsum <- sum(NewRecal_2$Ejweight) / 2^(nrow(NewRecal_2) - 1)
        NewRecal_2$masssum <- sum(NewRecal_2$mzweight) / 2^(nrow(NewRecal_2) - 1)
        Ejsum <- mean(NewRecal_2$Ejsum)
        peaks_out <- peaks_2
        peaks_out$mass <- peaks_out$mass / (1 + Ejsum)
        peaks_out <- rbind(peaks_final, peaks_out)
        Ej_out <- data.frame(Ejsum = Ejsum, Range = mean(NewRecal_2$Range))
        Ej_out <- rbind(Ej_final, Ej_out)
      }
    }
    peaks_final <- peaks_out
    Ej_final <- Ej_out
  }


  peaks <- peaks_final[peaks_final$mass > 0, ]
  Ej_final <- Ej_final[Ej_final$Range > 0, ]

  Ej <- Ej_final[Ej_final$Range == max(Ej_final$Range), ]

  peakshigh$mass <- peakshigh$mass / (1 + Ej$Ejsum)

  peaks <- rbind(peaks, peakshigh)
  # Separating the isotope and monoisotope masses for output.
  isopeaks2 <- subset(peaks, Tag2 != "Mono", select = c(1, 2, 3, 4))
  names(isopeaks2) <- c("Iso_Abund", "Iso_mass", "Iso_RT", "Tag2")

  peaks <- subset(peaks, Tag2 == "Mono", select = c(1, 2, 3))
  names(peaks) <- c("Abundance", "mass", "RT")

  ################################################
  # Recalibration for error comparison, calculating the average error before and after recalibration.

  if (cols == 3) {
    FinalRecal3 <- merge(df, FinalRecal2, by.x = c("Exp_mass", "Abundance", "RT", "formula"), by.y = c("Exp_mass", "Abundance", "RT", "formula"), all = T)
  }

  if (cols == 2) {
    FinalRecal3 <- merge(df, FinalRecal2, by.x = c("Exp_mass", "Abundance"), by.y = c("Exp_mass", "Abundance"), all = T)
  }
  FinalRecal3$Range <- 1
  for (i in 1:Dummy2) {
    I <- i

    FinalRecal3$Range <- replace(FinalRecal3$Range, FinalRecal3$Exp_mass >= min(FinalRecal3$Exp_mass) + (I - 1) *
      mzRange & FinalRecal3$Exp_mass < min(FinalRecal3$Exp_mass) + I * mzRange, I)
  }

  FinalRecal3$Range <- replace(FinalRecal3$Range, FinalRecal3$Range > Dummy3, Dummy3)

  FinalRecal3 <- merge(FinalRecal3, Ej_final, by.x = "Range", by.y = "Range")
  FinalRecal3$R_mass <- FinalRecal3$Exp_mass / (1 + FinalRecal3$Ejsum)
  FinalRecal3$R_err <- abs((FinalRecal3$R_mass - FinalRecal3$theor_mass) / FinalRecal3$theor_mass * 10^6)

  if (is.na(mean(FinalRecal3$Ejsum))) stop("Too many recalibrants in section,
                                         try decreasing CalPeak parameter") # New 5/13/20


  print(mean(FinalRecal3$AE_ppm))
  print(mean(FinalRecal3$R_err))


  plotpeak <- peaks[peaks$Abundance > SN, ]

  RecalOut <- merge(Align, Ej_final, by.x = "Range", by.y = "Range")
  RecalOut$R_mass <- RecalOut$Exp_mass / (1 + RecalOut$Ejsum)
  RecalOut$recal_err <- abs((RecalOut$R_mass - RecalOut$theor_mass) / RecalOut$theor_mass * 10^6)
  RecalOut$orig_err <- abs((RecalOut$Exp_mass - RecalOut$theor_mass) / RecalOut$theor_mass * 10^6)

  if (cols == 3) {
    AbundM <- df[c("Exp_mass", "Abundance", "RT")]
    AbundI <- df[c("C13_mass", "C13_abund", "C13_RT")]
    colnames(AbundI) <- c("Exp_mass", "Abundance", "RT")
    Abund <- rbind(AbundM, AbundI)
    Abund <- Abund[Abund$Abundance > 0, ]

    RecalPlot <- merge(RecalOut, Abund, by.x = c("Exp_mass", "Abundance", "RT"), by.y = c("Exp_mass", "Abundance", "RT"))

    RecalOut <- RecalOut[c(3, 5, 6, 9, 10)]

    # Preparing the final output of the function, the recalibrants list.
    peaks <- setNames(peaks[c(2, 1, 3)], c("exp_mass", "abundance", names(peaks)[3]))
    isopeaks2 <- setNames(isopeaks2[c(2, 1, 3, 4)], c("iso_mass", "iso_abund", names(isopeaks2)[3], "tag"))
    names(RecalOut)[c(1, 2)] <- c("abundance", "exp_mass")
  }

  if (cols == 2) {
    AbundM <- df[c("Exp_mass", "Abundance")]
    AbundI <- df[c("C13_mass", "C13_abund")]
    names(AbundI) <- c("Exp_mass", "Abundance")
    Abund <- rbind(AbundM, AbundI)
    Abund <- Abund[Abund$Abundance > 0, ]

    RecalPlot <- merge(RecalOut, Abund, by.x = c("Exp_mass", "Abundance"), by.y = c("Exp_mass", "Abundance"))

    RecalOut <- RecalOut[c(2, 3, 4, 7, 8)]

    # Preparing the final output of the function, the recalibrants list.
    peaks <- setNames(peaks[c(2, 1)], c("exp_mass", "abundance"))
    isopeaks2 <- setNames(isopeaks2[c(2, 1, 3)], c("iso_mass", "iso_abund", "tag"))
    names(RecalOut)[c(1, 2)] <- c("abundance", "exp_mass")
  }

  # Plot highlighting the recalibrant ions for qualitative assessment.
  MZ <- plot_spectrum(plotpeak, RecalPlot)

  Output <- list(Plot = MZ, Mono = peaks, Iso = isopeaks2, RecalList = RecalOut)
}

#' Identifies canditate series for recalibration
#'
#' RecalList() takes the output from MFAssign() and/or MFAssignCHO()
#' and identifies the homologous series that could be used as recalibrants.
#'
#' It returns a dataframe that contains the CH2 homologous series that
#' contain more than 3 members.
#'
#' The columns of the returned dataframe are as follows:
#'
#' Series - reports the homologous series according to class, adduct, and DBE,
#' the format is "class_adduct_DBE", for example a homologous series with
#' class = "O6, adduct of Na+, and DBE = 4 would be "O6_Na_4"
#'
#' Number Observed - reports the number of members of each homologous series.
#'
#' Series Index - represents the order of the series when ordered by length of
#' homologous series.
#'
#' Mass Range - reports the minimum and maximum mass for the compounds within a
#' homologous series.
#'
#' Tall Peak - reports the mass of the most abundant peak in each series.
# ’
#' Abundance Score - reports the percentage difference between the mean abundance
#' of a homologous series and the median abundance within the mass range the
#' "Tall Peak" falls in (for example m/z 200-300). A higher score is generally better.
#'
#' Peak Score - This column compares the intensity of the tallest peak in a
#' given series to the second tallest peak in the series This comparison is
#' calculated by log10(Max Peak Intensity/Second Peak Intensity) The closer
#' to 0 this value is the better, in general.
#'
#' Peak Distance - This column shows the number of CH2 units between the tallest
#' and second tallest peak in each series. In general, it is better for the
#' value to be as close to 1 as possible.
#'
#' Series Score - This column compares the number of actual observations in each
#' series to the theoretical maximum number based on the CH2 homologous series.
#' The closer to one this value is, the better.
#'
#'
#' While the function does make some minor decisions as to which series to choose
#' from (more than 3 members of homologous series) the rest of the decision
#' making is left to the user.
#'
#' @param df data frame output from MFAssign() or MFAssignCHO()
#' @return data frame
#' @export
RecalList <- function(df) {
  cols_to_remove_if_49 <- c(19:21)
  cols_to_remove_if_53 <- c(3, 20:22, 46, 49, 52)

  df <- df[, -which(
    ncol(df) == 49 & seq_along(df) %in% cols_to_remove_if_49 |
      ncol(df) == 53 & seq_along(df) %in% cols_to_remove_if_53
  ), drop = FALSE]


  df$number <- 1
  df$Adduct <- ifelse(df$M > 0, "Na", "H")
  df$Adduct[df$POE == 1 | df$NOE == 1] <- "OE"
  df$SeriesAdd <- paste(df$class, df$Adduct, sep = "_")

  colnames(df) <- sub("exp_mass", "Exp_mass", colnames(df))
  colnames(df) <- sub("abundance", "Abundance", colnames(df))

  df1 <- na.omit(aggregate(number ~ SeriesAdd + DBE, df, sum, na.rm = TRUE))

  longseries <- df1 |>
    dplyr::arrange(desc(number)) |>
    dplyr::mutate(Index = dplyr::row_number())
  Recal <- merge(df, longseries, by.x = c("SeriesAdd", "DBE"), by.y = c("SeriesAdd", "DBE"))
  Recal <- tidyr::unite(Recal, Series, SeriesAdd, DBE, sep = "_", remove = FALSE)

  Recal <- AddCalculatedSummary(Recal)
  Recal <- FilterAndMerge(Recal, c(56), "Maxmass")
  Recal <- FilterAndMerge(Recal, c(59), "Secmass")

  Recal <- Recal[Recal$number.y > 3, ]
  names(Recal)[56] <- "Maxmass"
  Recal <- Recal[!is.na(Recal$Maxmass), ]

  breaks <- c(0, 200, 300, 400, 500, 600, 700, 800, Inf)
  Recal$M.window <- cut(Recal$Maxmass, breaks,
    labels = c("0-200", "200-300", "300-400", "400-500", "500-600", "600-700", "700-800", ">800"),
    include.lowest = TRUE
  )

  Recal1 <- aggregate(MInt ~ M.window, Recal, median, na.rm = TRUE)
  Recal <- merge(Recal, Recal1, by = "M.window")

  Recal$IntRel <- (Recal$MInt.x - Recal$MInt.y) / Recal$MInt.y * 100
  Recal$Peakcomp <- log10(Recal$Maxint / Recal$Secint)
  Recal$Nextpeak <- abs((Recal$Maxmass.y - Recal$Secmass.y) / 14)
  Recal <- Recal[c(2, 52:66)]
  Recal$Min <- round(Recal$Min, 3)
  Recal$Max <- round(Recal$Max, 3)
  Recal$SerScor <- ((Recal$Max - Recal$Min) / 14 + 1) / Recal$number.y
  Recal <- tidyr::unite(Recal, Range, Min, Max, sep = "-")

  columns_to_drop <- c(5, 7, 8, 9, 10, 11, 12)
  Recal <- Recal[, -columns_to_drop]

  columns_to_rename <- c(2, 3, 4, 5, 6, 7, 8, 9)
  new_column_names <- c("Number Observed", "Series Index", "Mass Range", "Tall Peak", "Abundance Score", "Peak Score", "Peak Distance", "Series Score")
  names(Recal)[columns_to_rename] <- new_column_names
  Recal <- Recal[!duplicated(Recal), ]
}

#' AddCalculatedSummary function
#'
#' This function takes a data frame 'Recal' and calculates summary statistics
#' grouped by the 'Index' column. It adds columns for Minimum, Maximum, Mean,
#' Maximum Mass, Maximum Intensity, Second Maximum Intensity, and Second Maximum Mass.
#'
#' @param Recal data frame
#' @return A modified data frame with added calculated summary columns.
#'
AddCalculatedSummary <- function(Recal) {
  Recal <- dplyr::group_by(Recal, Index)
  Recal <- dplyr::mutate(Recal,
    Min = min(Exp_mass),
    Max = max(Exp_mass),
    MInt = mean(Abundance),
    Maxmass = ifelse(Abundance == max(Abundance),
      Exp_mass, NA
    ),
    Maxint = max(Abundance),
    Secint = sort(Abundance, TRUE)[2],
    Secmass = ifelse(Abundance == sort(Abundance, TRUE)[2], Exp_mass, NA)
  )
  return(Recal)
}


#' FilterAndMerge function
#'
#' This function takes a data frame 'Recal', a vector of column indices 'column_indices',
#' and a column name 'column_name'. It performs the following steps:
#' 1. Selects specified columns from 'Recal'.
#' 2. Filters rows with complete cases for the specified column.
#' 3. Merges the original 'Recal' data frame with the filtered columns.
#'
#' @param Recal The input data frame.
#' @param column_indices A vector of column indices to be selected.
#' @param column_name The name of the column used for filtering.
#'
#' @return  A modified data frame with selected columns and filtered rows.
#'
FilterAndMerge <- function(Recal, column_indices, column_name) {
  selected_columns <- Recal[, c(1, column_indices)]
  selected_columns <- selected_columns[complete.cases(selected_columns[, column_name]), ]
  Recal <- merge(Recal, selected_columns, by.x = c("Series"), by.y = c("Series"))
  return(Recal)
}
