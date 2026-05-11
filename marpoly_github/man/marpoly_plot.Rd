\name{marpoly_plot}
\alias{marpoly_plot}
\title{Quick MAR and polyspectral plots}
\usage{marpoly_plot(x, orders = c(2, 3), fit = NULL, grid_size = NULL,
  component = c("real", "imaginary", "modulus"), slice = 1L, demean = TRUE)}
\arguments{
\item{x}{Numeric time series.}
\item{orders}{Spectral orders to display.}
\item{fit}{Optional fitted \code{marpoly_fit} object.}
\item{grid_size}{Frequency grid size.}
\item{component}{Component for complex heatmaps.}
\item{slice}{Order-4 slice.}
\item{demean}{Whether to demean.}
}
\description{Plots a series, selected periodograms, and optional fitted spectrum.}
\value{Invisibly returns periodogram objects.}
