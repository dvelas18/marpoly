\name{marpoly_simulate}
\alias{marpoly_simulate}
\title{Simulate a MAR(r,s) process}
\usage{marpoly_simulate(n, causal = NULL, noncausal = NULL, innov = NULL,
  dist = "normal", burnin = 500, seed = NULL, standardize = TRUE,
  check_stationary = TRUE, ...)}
\arguments{
\item{n}{Number of observations.}
\item{causal}{Causal lag coefficients.}
\item{noncausal}{Noncausal lead coefficients.}
\item{innov}{Optional vector or function generating innovations.}
\item{dist}{Innovation distribution: normal, t, skewt, or a function.}
\item{burnin}{Number of observations discarded from each side.}
\item{seed}{Optional seed.}
\item{standardize}{Whether to center and scale innovations.}
\item{check_stationary}{Whether to enforce roots outside the unit circle.}
\item{...}{Additional innovation arguments such as df and gamma.}
}
\description{Simulates \code{(1 - causal(L))(1 - noncausal(L^-1)) y_t = epsilon_t}.}
\value{A \code{marpoly_sim} object.}
