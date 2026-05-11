\name{marpoly_diagnostics}
\alias{marpoly_diagnostics}
\alias{MAR_error_diagnostics}
\title{Error diagnostics}
\usage{
marpoly_diagnostics(errors, powers = c(1L, 2L), lag.max = 20L, lags = NULL)
MAR_error_diagnostics(errors, lags = c(1L, 2L, 5L, 10L), powers = 1:4)
}
\arguments{
\item{errors}{Numeric fitted-error vector.}
\item{powers}{Powers of errors to diagnose.}
\item{lag.max}{Maximum lag included by the wrapper.}
\item{lags}{Ljung-Box lags.}
}
\description{Computes summary statistics and Ljung-Box tests for fitted errors and powers of fitted errors.}
\value{A list with a summary table and a Ljung-Box table.}
