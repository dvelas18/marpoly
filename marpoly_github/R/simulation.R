#' Simulate a MAR(r,s) process
#'
#' Simulates the mixed causal-noncausal autoregressive model
#' `(1 - causal(L))(1 - noncausal(L^-1)) y_t = epsilon_t`.
#' The simulator first solves the causal component forward in time and then
#' solves the noncausal component by applying the same recursion to the reversed
#' series.
#'
#' @param n Number of observations to return.
#' @param causal Numeric vector of causal lag coefficients.
#' @param noncausal Numeric vector of noncausal lead coefficients.
#' @param innov Optional numeric vector or function generating innovations. If a
#'   function, it must accept one argument, the requested sample size.
#' @param dist Innovation distribution used when `innov` is not supplied.
#' @param burnin Number of observations discarded from each side.
#' @param seed Optional random seed.
#' @param standardize If TRUE, center and scale generated innovations.
#' @param check_stationary If TRUE, enforce the MAR stationarity restrictions.
#' @param ... Additional arguments passed to `marpoly_rinnov`, such as `df` and
#'   `gamma`.
#' @return An object of class `marpoly_sim` containing the simulated series,
#'   innovations, and model parameters.
#' @export
marpoly_simulate <- function(n, causal = NULL, noncausal = NULL, innov = NULL,
                             dist = "normal", burnin = 500, seed = NULL,
                             standardize = TRUE, check_stationary = TRUE, ...) {
  n <- as.integer(n)
  burnin <- as.integer(burnin)
  if (n <= 0L) stop("n must be positive.", call. = FALSE)
  if (burnin < 0L) stop("burnin must be nonnegative.", call. = FALSE)
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  if (check_stationary) {
    mp_check_stationary(causal, noncausal)
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
  }
  total <- n + 2L * burnin
  if (is.null(innov)) {
    eps <- marpoly_rinnov(total, dist = dist, standardize = standardize, ...)
  } else if (is.function(innov)) {
    eps <- innov(total)
    eps <- as.numeric(eps)
    if (standardize) {
      sx <- stats::sd(eps)
      if (is.finite(sx) && sx > 0) eps <- (eps - mean(eps)) / sx else eps <- eps - mean(eps)
    }
  } else {
    eps <- as.numeric(innov)
    if (length(eps) < total) {
      stop("innov must have at least n + 2 * burnin observations.", call. = FALSE)
    }
    eps <- eps[seq_len(total)]
  }
  if (any(!is.finite(eps))) stop("Innovations must be finite.", call. = FALSE)

  causal_part <- mp_ar_filter_forward(eps, causal)
  mixed_reversed <- mp_ar_filter_forward(rev(causal_part), noncausal)
  y_full <- rev(mixed_reversed)
  keep <- (burnin + 1L):(burnin + n)
  y <- y_full[keep]
  eps_keep <- eps[keep]
  out <- list(
    x = y,
    innovations = eps_keep,
    full_x = y_full,
    full_innovations = eps,
    causal = causal,
    noncausal = noncausal,
    n = n,
    burnin = burnin,
    model = mp_model_label(length(causal), length(noncausal)),
    stationary = marpoly_stationary(causal, noncausal)
  )
  class(out) <- "marpoly_sim"
  out
}

#' @export
print.marpoly_sim <- function(x, ...) {
  cat("marpoly simulation\n")
  cat("  model:      ", x$model, "\n", sep = "")
  cat("  n:          ", x$n, "\n", sep = "")
  cat("  causal:     ", if (length(x$causal)) mp_format_vector(x$causal) else "none", "\n", sep = "")
  cat("  noncausal:  ", if (length(x$noncausal)) mp_format_vector(x$noncausal) else "none", "\n", sep = "")
  invisible(x)
}
