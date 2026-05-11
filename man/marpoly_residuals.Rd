\name{marpoly_residuals}
\alias{marpoly_residuals}
\alias{marpoly_errors}
\title{MAR residuals and fitted errors}
\usage{marpoly_residuals(x, causal = NULL, noncausal = NULL, drop = TRUE, demean = TRUE)
marpoly_errors(fit, allocation = "selected", demean = TRUE)}
\arguments{
\item{x}{Numeric time series.}
\item{causal}{Causal coefficients.}
\item{noncausal}{Noncausal coefficients.}
\item{drop}{Whether to drop boundary observations.}
\item{demean}{Whether to demean before filtering.}
\item{fit}{A \code{marpoly_fit} object.}
\item{allocation}{selected, all, or an allocation id.}
}
\description{Computes MAR errors for supplied coefficients or extracts errors from fitted allocations.}
\value{A residual vector or a list of residual vectors.}
