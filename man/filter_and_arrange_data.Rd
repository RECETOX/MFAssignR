% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IsoFiltR.R
\name{filter_and_arrange_data}
\alias{filter_and_arrange_data}
\title{Filters and arranges data based on abundance and mass}
\usage{
filter_and_arrange_data(peaks, sn)
}
\arguments{
\item{peaks}{data frame:
The input data frame containing at least two columns, where the first column is the experimental mass,
and the second column is the abundance.}

\item{sn}{numeric:
The abundance threshold. Peaks with abundance below this threshold will be filtered out.}
}
\value{
data frame:
A filtered and arranged data frame with columns "Exp_mass" and "Abundance".
}
\description{
This function renames columns and filters data based on a specified abundance threshold.
It arranges the data based on the experimental mass.
}
