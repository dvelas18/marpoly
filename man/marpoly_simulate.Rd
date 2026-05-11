\name{marpoly_simulate}
\alias{marpoly_simulate}
\alias{simulate_MAR}
\alias{rskt}
\title{Simulate a MAR(r,s) process and skew Student t innovations}
\usage{
marpoly_simulate(n = NULL, causal = NULL, noncausal = NULL, innov = NULL,
  dist = "normal", burnin = 500, seed = NULL, standardize = TRUE,
  check_stationary = TRUE, ..., T = NULL, phi_plus = NULL,
  phi_star = NULL, noise_dist = NULL, noise_args = list())
simulate_MAR(T, phi_plus = NULL, phi_star = NULL, noise_dist = stats::rnorm,
  noise_args = list(), burnin = 500, seed = NULL, standardize = TRUE,
  check_stationarity = TRUE)
rskt(n, df = 5, gamma = 1, standardize = TRUE)
}
\arguments{
\item{n, T}{Number of observations.}
\item{causal, phi_plus}{Causal lag coefficients.}
\item{noncausal, phi_star}{Noncausal lead coefficients.}
\item{innov}{Optional vector or function generating innovations.}
\item{dist}{Innovation distribution: \code{"normal"}, \code{"t"}, or \code{"skewt"}.}
\item{noise_dist}{Innovation generator function for \code{simulate_MAR}.}
\item{noise_args}{List of arguments passed to \code{noise_dist}.}
\item{burnin}{Number of observations discarded from each side.}
\item{seed}{Optional seed.}
\item{standardize}{Whether to center and scale innovations.}
\item{check_stationary, check_stationarity}{Whether to enforce roots outside the unit circle.}
\item{df}{Degrees of freedom for Student t or skew Student t innovations.}
\item{gamma}{Fernandez-Steel skewness parameter.}
\item{...}{Additional innovation arguments such as \code{df} and \code{gamma}.}
}
\description{Simulates \code{(1 - phi_plus(L))(1 - phi_star(L^-1)) y_t = epsilon_t}. Both component polynomials are checked for roots outside the unit circle unless stationarity checking is disabled.}
\value{A \code{marpoly_sim} object. The series is available as both \code{x} and \code{X}.}
