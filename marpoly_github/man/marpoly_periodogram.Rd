\name{marpoly_periodogram}
\alias{marpoly_periodogram}
\alias{marpoly_dft}
\alias{marpoly_periodograms}
\title{Higher-order periodograms}
\usage{marpoly_dft(x, demean = TRUE)
marpoly_periodogram(x, order = 2, grid_size = NULL, demean = TRUE, indices = NULL)
marpoly_periodograms(x, orders = c(2, 3, 4), ...)}
\arguments{
\item{x}{Numeric time series.}
\item{order}{Order 2, 3, or 4.}
\item{grid_size}{Number of frequencies in each dimension.}
\item{demean}{Whether to demean the series.}
\item{indices}{Optional Fourier indices.}
\item{orders}{Vector of orders.}
\item{...}{Arguments passed to \code{marpoly_periodogram}.}
}
\description{Computes the discrete Fourier transform, periodogram, biperiodogram, and a fourth-order periodogram.}
\value{A \code{marpoly_periodogram} object or a list of such objects.}
