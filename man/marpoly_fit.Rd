\name{marpoly_fit}
\alias{marpoly_fit}
\alias{marpoly_fit_all}
\alias{marpoly_rank}
\alias{marpoly_order_sets}
\title{Estimate and rank MAR root allocations}
\usage{marpoly_fit(x, p = NULL, p_max = 6, orders = c(2, 3), weights = NULL,
  grid_size = NULL, allocations = "all", optimizer = c("single", "sann"),
  maxit = 3000, criterion = c("BIC", "AIC"), demean = TRUE,
  include_mean = FALSE, seed = NULL, se = TRUE)
marpoly_fit_all(...)
marpoly_rank(fit, n = 10)
marpoly_order_sets(orders = c(2, 3, 4))}
\arguments{
\item{x}{Numeric time series.}
\item{p}{Pseudo-causal AR order. If NULL, selected by BIC/AIC.}
\item{p_max}{Maximum order when selecting p.}
\item{orders}{Spectral orders used in the objective.}
\item{weights}{Objective weights.}
\item{grid_size}{Frequency grid size.}
\item{allocations}{Root allocations to estimate.}
\item{optimizer}{Local single-start optimization or simulated annealing.}
\item{maxit}{Maximum optimizer iterations.}
\item{criterion}{Information criterion for p selection.}
\item{demean}{Whether to demean the input series.}
\item{include_mean}{Whether to include intercept in preliminary AR regression.}
\item{seed}{Optional seed.}
\item{se}{Whether to compute numerical-Hessian standard errors.}
\item{...}{Arguments passed to \code{marpoly_fit}.}
\item{fit}{A \code{marpoly_fit} object.}
\item{n}{Rows to return in rank table.}
}
\description{Estimates all selected MAR root allocations and ranks them by a frequency-domain minimum-distance objective based on spectral orders 2, 3, and/or 4.}
\value{A \code{marpoly_fit} object or a rank table.}
