\name{marpoly_fit}
\alias{marpoly_fit}
\alias{marpoly_fit_all}
\alias{marpoly_rank}
\alias{marpoly_order_sets}
\title{Estimate and rank MAR models by polyspectral minimum distance}
\usage{
marpoly_fit(x, p = NULL, p_max = 6, orders = c(2, 3), weights = NULL,
  grid_size = NULL, allocations = "all", optimizer = c("single", "sann"),
  maxit = 3000, criterion = c("BIC", "AIC", "AICc"), demean = TRUE,
  include_mean = FALSE, seed = NULL, se = TRUE, control = list(),
  ar_method = c("yule-walker", "burg", "ols", "mle"), ...)
marpoly_fit_all(...)
marpoly_rank(fit, n = 10, root_allocations = FALSE)
marpoly_order_sets(orders = c(2, 3, 4))
}
\arguments{
\item{x}{Numeric time series.}
\item{p}{Fixed pseudo-causal AR order. If \code{NULL}, it is selected by second-order information criterion.}
\item{p_max}{Maximum AR order when \code{p} is selected automatically.}
\item{orders}{Spectral orders used in the objective: any subset of \code{c(2,3,4)}.}
\item{weights}{Optional objective weights, in the same order as \code{orders}.}
\item{grid_size}{Compatibility shortcut for coarser fourth-order frequency selection. For full control, use \code{control}.}
\item{allocations}{MAR(r,s) specifications to estimate when \code{p} is fixed. Use \code{"all"}, model labels, a two-column data frame, or a list of \code{c(r,s)} pairs.}
\item{optimizer}{\code{"single"} for root-allocation starts with Nelder-Mead, or \code{"sann"} for simulated annealing.}
\item{maxit}{Optimizer iteration limit.}
\item{criterion}{Information criterion for automatic p selection.}
\item{demean}{If TRUE, demean the series.}
\item{include_mean}{Kept for backward compatibility; the spectral workflow controls centering through \code{demean}.}
\item{seed}{Optional random seed.}
\item{se}{If TRUE, compute finite-difference sandwich standard errors.}
\item{control}{Advanced control list passed to \code{MAR_spectral_estimate}.}
\item{ar_method}{Preliminary causal AR method: \code{"yule-walker"}, \code{"burg"}, \code{"ols"}, or \code{"mle"}.}
\item{...}{Additional control arguments.}
\item{fit}{A fitted object returned by \code{marpoly_fit}.}
\item{n}{Rows to return in rank table.}
\item{root_allocations}{If TRUE, rank every root allocation/local solution instead of one row per MAR(r,s).}
}
\description{
Estimates MAR(r,s) models with the spectral, bispectral, and trispectral minimum-distance objective.
The implementation first estimates a pseudo-causal AR(p), constructs root-allocation starts, and then ranks the feasible causal/noncausal specifications by the selected spectral objective.
}
\value{A \code{marpoly_fit} object, or a rank table for \code{marpoly_rank}.}
