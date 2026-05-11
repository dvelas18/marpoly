\name{marpoly_spectrum}
\alias{marpoly_spectrum}
\alias{marpoly_spectra}
\alias{marpoly_transfer}
\title{MAR transfer function and polyspectral densities}
\usage{marpoly_transfer(omega, causal = NULL, noncausal = NULL)
marpoly_spectrum(causal = NULL, noncausal = NULL, order = 2, kappa = 1,
  grid_size = NULL, omega = NULL, n = 512)
marpoly_spectra(causal = NULL, noncausal = NULL, orders = c(2, 3, 4), kappa = NULL, ...)}
\arguments{
\item{omega}{Angular frequencies.}
\item{causal}{Causal coefficients.}
\item{noncausal}{Noncausal coefficients.}
\item{order}{Order 2, 3, or 4.}
\item{kappa}{Innovation cumulant.}
\item{grid_size}{Number of frequencies in each dimension.}
\item{n}{Frequency denominator when building a grid.}
\item{orders}{Vector of orders.}
\item{...}{Additional arguments.}
}
\description{Computes MAR transfer functions and theoretical spectrum, bispectrum, or fourth-order spectral density.}
\value{A complex vector/array or a \code{marpoly_spectrum} object.}
