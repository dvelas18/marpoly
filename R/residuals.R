mp_lag_coefficients <- function(causal = NULL, noncausal = NULL) {
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  r <- length(causal)
  s <- length(noncausal)
  lags <- -s:r
  vals <- setNames(rep(0, length(lags)), as.character(lags))
  vals[as.character(0)] <- vals[as.character(0)] + 1
  if (r > 0L) {
    for (i in seq_len(r)) vals[as.character(i)] <- vals[as.character(i)] - causal[i]
  }
  if (s > 0L) {
    for (j in seq_len(s)) vals[as.character(-j)] <- vals[as.character(-j)] - noncausal[j]
  }
  if (r > 0L && s > 0L) {
    for (i in seq_len(r)) {
      for (j in seq_len(s)) {
        lag <- i - j
        vals[as.character(lag)] <- vals[as.character(lag)] + causal[i] * noncausal[j]
      }
    }
  }
  data.frame(lag = as.integer(names(vals)), coefficient = as.numeric(vals), row.names = NULL)
}

#' Extract MAR residuals for supplied coefficients
#'
#' @param x Numeric time series.
#' @param causal Numeric vector of causal lag coefficients.
#' @param noncausal Numeric vector of noncausal lead coefficients.
#' @param drop If TRUE, return only observations for which all lags/leads are
#'   observed. If FALSE, use NA at the boundaries.
#' @param demean If TRUE, subtract the sample mean before filtering.
#' @return Numeric vector of residuals. When `drop = TRUE`, the returned vector
#'   has attribute `time_index` giving the original time indices.
#' @export
marpoly_residuals <- function(x, causal = NULL, noncausal = NULL, drop = TRUE,
                              demean = TRUE) {
  x <- mp_clean_series(x, demean = demean)
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  n <- length(x)
  r <- length(causal)
  s <- length(noncausal)
  lag_table <- mp_lag_coefficients(causal, noncausal)
  eps <- rep(NA_real_, n)
  start <- if (r > 0L) r + 1L else 1L
  end <- if (s > 0L) n - s else n
  if (start > end) {
    stop("Not enough observations to compute residuals for this MAR order.", call. = FALSE)
  }
  for (t in start:end) {
    val <- 0
    for (j in seq_len(nrow(lag_table))) {
      lag <- lag_table$lag[j]
      val <- val + lag_table$coefficient[j] * x[t - lag]
    }
    eps[t] <- val
  }
  if (drop) {
    out <- eps[start:end]
    attr(out, "time_index") <- start:end
    return(out)
  }
  eps
}

#' Extract residuals/errors from a fitted MAR allocation
#'
#' @param object A `marpoly_fit` object.
#' @param allocation "selected", "all", or an allocation id.
#' @param ... Reserved for S3 compatibility.
#' @return Numeric vector for one allocation or a list of vectors for all
#'   allocations.
#' @export
residuals.marpoly_fit <- function(object, allocation = "selected", ...) {
  marpoly_errors(object, allocation = allocation, ...)
}

#' Extract fitted MAR errors for selected or all allocations
#'
#' @param fit A `marpoly_fit` object.
#' @param allocation "selected", "all", or an allocation id.
#' @param demean If TRUE, demean before filtering.
#' @return Residual vector or named list of residual vectors.
#' @export
marpoly_errors <- function(fit, allocation = "selected", demean = TRUE) {
  if (!inherits(fit, "marpoly_fit")) {
    stop("fit must be a marpoly_fit object.", call. = FALSE)
  }
  if (identical(allocation, "selected")) {
    row <- fit$results[fit$selected_index, ]
    return(marpoly_residuals(fit$x, causal = row$causal[[1L]], noncausal = row$noncausal[[1L]], demean = demean))
  }
  if (identical(allocation, "all")) {
    out <- vector("list", nrow(fit$results))
    for (i in seq_len(nrow(fit$results))) {
      out[[i]] <- marpoly_residuals(fit$x, causal = fit$results$causal[[i]], noncausal = fit$results$noncausal[[i]], demean = demean)
    }
    names(out) <- paste0(fit$results$rank, "_", fit$results$model, "_allocation", fit$results$allocation_id)
    return(out)
  }
  id <- as.integer(allocation)
  pos <- which(fit$results$allocation_id == id)
  if (length(pos) != 1L) stop("allocation id not found.", call. = FALSE)
  marpoly_residuals(fit$x, causal = fit$results$causal[[pos]], noncausal = fit$results$noncausal[[pos]], demean = demean)
}
