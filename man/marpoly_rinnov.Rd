\name{marpoly_rinnov}
\alias{marpoly_rinnov}
\alias{marpoly_rskewt}
\title{Innovation generators}
\usage{marpoly_rinnov(n, dist = "normal", df = 5, gamma = 1, standardize = TRUE, ...)
marpoly_rskewt(n, df = 5, gamma = 1, standardize = TRUE)}
\arguments{
\item{n}{Number of observations.}
\item{dist}{Distribution name or generator function.}
\item{df}{Degrees of freedom.}
\item{gamma}{Skewness parameter.}
\item{standardize}{Whether to center and scale the sample.}
\item{...}{Additional arguments.}
}
\description{Generates normal, Student t, or Fernandez-Steel skew Student t innovations.}
\value{Numeric vector.}
