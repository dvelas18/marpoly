\name{MAR_spectral_estimate}
\alias{MAR_spectral_estimate}
\alias{MAR_select_all_by_second_order}
\alias{MAR_select_all_fixed_p}
\alias{MAR_select_all_by_auto_arima}
\alias{MAR_select_by_order_sets}
\alias{MAR_all_order_sets}
\alias{MAR_spectrum_indices}
\alias{MAR_bispectrum_pairs}
\alias{MAR_trispectrum_tuples}
\alias{MAR_estimate_causal_eta}
\alias{MAR_select_p_second_order}
\alias{MAR_all_specs_from_p}
\alias{MAR_auto_starts_from_eta}
\alias{MAR_model_errors}
\alias{MAR_objective_prepare}
\alias{MAR_standard_errors}
\alias{MAR_root_combination_table}
\alias{MAR_best_models_table}
\alias{MAR_report}
\title{MAR spectral minimum-distance estimation workflow}
\usage{
MAR_spectral_estimate(y, r, s, orders = c(2, 3), order = NULL,
  start = NULL, control = list(), ...)
MAR_select_all_by_second_order(y, orders = c(2, 3), order = NULL,
  max_p = 8L, ic = c("bic", "aic", "aicc"),
  ar_method = c("yule-walker", "burg", "ols", "mle"),
  specs = NULL, control = list(), ...)
MAR_select_all_fixed_p(y, p, orders = c(2, 3),
  ar_method = c("yule-walker", "burg", "ols", "mle"),
  specs = NULL, control = list(), allocations = "all", ...)
simulate_MAR(T, phi_plus = NULL, phi_star = NULL, noise_dist = NULL,
  noise_args = list(), seed = NULL, burnin = 500, standardize = TRUE,
  check_stationary = TRUE, ...)
rskt(n, df = 5, gamma = 1, standardize = TRUE)
}
\arguments{
\item{y}{Numeric time series.}
\item{r}{Causal MAR order.}
\item{s}{Noncausal MAR order.}
\item{orders, order}{Spectral orders to include: 2, 3, 4, or any subset.}
\item{start}{Optional starting values.}
\item{control}{Control list for frequency selection, optimization, errors, and standard errors.}
\item{p}{Fixed pseudo-causal AR order.}
\item{max_p}{Maximum preliminary pseudo-causal AR order.}
\item{ic}{Information criterion.}
\item{ar_method}{Preliminary AR estimation method passed to \code{stats::ar}.}
\item{specs}{Optional data frame with columns \code{r} and \code{s}.}
\item{allocations}{Optional specification filter for the fixed-p workflow.}
\item{T}{Number of simulated observations.}
\item{phi_plus}{Causal coefficients.}
\item{phi_star}{Noncausal coefficients.}
\item{noise_dist}{Innovation generator.}
\item{noise_args}{Arguments passed to \code{noise_dist}.}
\item{seed}{Optional random seed.}
\item{burnin}{Burn-in observations discarded from both sides.}
\item{standardize}{Whether to center and scale generated innovations.}
\item{check_stationary}{Whether to enforce stationarity restrictions.}
\item{n}{Number of random draws.}
\item{df}{Degrees of freedom for the skew Student t generator.}
\item{gamma}{Skewness parameter for the skew Student t generator.}
\item{...}{Additional arguments.}
}
\description{
Functions implementing the MAR spectral, bispectral, and trispectral minimum-distance estimation workflow. The workflow selects a pseudo-causal AR order by second-order criteria, constructs root-allocation starting values, estimates MAR(r,s) specifications using selected spectral orders, ranks the selected model and competitors, and optionally computes fitted errors and diagnostics.
}
\value{A fitted MAR object, selection object, data frame, vector, or list depending on the function.}
