\name{marpoly_diagnostics}
\alias{marpoly_diagnostics}
\title{Residual diagnostics}
\usage{marpoly_diagnostics(residuals, lag.max = 20, powers = c(1, 2), fitdf = 0)}
\arguments{
\item{residuals}{Numeric residual vector.}
\item{lag.max}{Maximum lag.}
\item{powers}{Powers of residuals to diagnose.}
\item{fitdf}{Degrees of freedom in Ljung-Box tests.}
}
\description{Computes ACF values and Ljung-Box tests for residual powers.}
\value{A \code{marpoly_diagnostics} object.}
