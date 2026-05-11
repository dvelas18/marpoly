\name{marpoly_stationary}
\alias{marpoly_stationary}
\title{Check MAR stationarity}
\usage{marpoly_stationary(causal = NULL, noncausal = NULL, tol = 1e-8)}
\arguments{
\item{causal}{Causal coefficients.}
\item{noncausal}{Noncausal coefficients.}
\item{tol}{Unit-root tolerance.}
}
\description{Checks whether both MAR component polynomials have roots outside the unit circle.}
\value{A list with the stationarity result and roots.}
