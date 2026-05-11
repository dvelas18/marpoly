mp_prepare_periodograms <- function(x, orders, grid_size = NULL, demean = TRUE) {
  orders <- mp_parse_orders(orders)
  out <- list()
  for (k in orders) {
    gs <- if (is.null(grid_size)) NULL else if (length(grid_size) == 1L) grid_size else grid_size[as.character(k)]
    if (length(gs) == 0L || is.na(gs)) gs <- NULL
    out[[as.character(k)]] <- marpoly_periodogram(x, order = k, grid_size = gs, demean = demean)
  }
  out
}

mp_kappa_estimate <- function(Ik, prod_theta, order) {
  raw <- Re(Ik / prod_theta)
  raw <- raw[is.finite(raw)]
  if (length(raw) == 0L) return(0)
  (2 * pi)^(order - 1L) * mean(raw)
}

mp_objective_components <- function(causal, noncausal, pdata, eta, orders, weights) {
  total <- 0
  kappas <- numeric(length(orders))
  names(kappas) <- as.character(orders)
  for (kk in seq_along(orders)) {
    k <- orders[kk]
    pg <- pdata[[as.character(k)]]
    prod_theta <- mp_transfer_product(pg$omega, order = k, causal = causal, noncausal = noncausal)
    prod_eta <- mp_transfer_product(pg$omega, order = k, causal = eta, noncausal = NULL)
    kappa <- mp_kappa_estimate(pg$values, prod_theta, order = k)
    kappas[kk] <- kappa
    Sk <- kappa / (2 * pi)^(k - 1L) * prod_theta
    if (k == 2L) Sk <- Re(Sk)
    den <- Mod(prod_eta)^2
    den[!is.finite(den) | den < 1e-12] <- 1e-12
    part <- mean(Mod(pg$values - Sk)^2 / den, na.rm = TRUE)
    if (!is.finite(part)) part <- 1e12
    total <- total + weights[kk] * part
  }
  list(value = total, kappas = kappas)
}

mp_objective <- function(par, r, s, pdata, eta, orders, weights, penalty = 1e12) {
  causal <- if (r > 0L) par[seq_len(r)] else numeric(0)
  noncausal <- if (s > 0L) par[r + seq_len(s)] else numeric(0)
  if (any(!is.finite(par))) return(penalty)
  if (!mp_is_stationary_poly(causal, tol = -1e-5) || !mp_is_stationary_poly(noncausal, tol = -1e-5)) {
    return(penalty + 1e6 * sum(par^2))
  }
  out <- mp_objective_components(causal, noncausal, pdata, eta, orders, weights)$value
  if (!is.finite(out)) penalty else out
}

mp_fit_one_allocation <- function(x, allocation_row, pdata, eta, orders, weights,
                                  optimizer = "single", maxit = 3000, sann_temp = 10,
                                  polish = TRUE, se = TRUE) {
  r <- allocation_row$r
  s <- allocation_row$s
  start <- c(allocation_row$causal[[1L]], allocation_row$noncausal[[1L]])
  obj <- function(par) mp_objective(par, r = r, s = s, pdata = pdata, eta = eta, orders = orders, weights = weights)
  if (length(start) == 0L) {
    par_hat <- numeric(0)
    value <- obj(par_hat)
    convergence <- 0L
    message <- "white-noise"
  } else {
    method <- if (optimizer == "sann") "SANN" else "Nelder-Mead"
    control <- list(maxit = maxit)
    if (optimizer == "sann") control$temp <- sann_temp
    opt <- try(stats::optim(start, obj, method = method, control = control), silent = TRUE)
    if (inherits(opt, "try-error")) {
      par_hat <- start
      value <- obj(start)
      convergence <- 999L
      message <- as.character(opt)
    } else {
      par_hat <- opt$par
      value <- opt$value
      convergence <- opt$convergence
      message <- opt$message
    }
    if (optimizer == "sann" && polish && length(par_hat) > 0L) {
      opt2 <- try(stats::optim(par_hat, obj, method = "Nelder-Mead", control = list(maxit = max(500L, floor(maxit / 3L)))), silent = TRUE)
      if (!inherits(opt2, "try-error") && is.finite(opt2$value) && opt2$value <= value) {
        par_hat <- opt2$par
        value <- opt2$value
        convergence <- opt2$convergence
        message <- opt2$message
      }
    }
  }
  causal_hat <- if (r > 0L) par_hat[seq_len(r)] else numeric(0)
  noncausal_hat <- if (s > 0L) par_hat[r + seq_len(s)] else numeric(0)
  comps <- mp_objective_components(causal_hat, noncausal_hat, pdata, eta, orders, weights)
  se_vec <- rep(NA_real_, length(par_hat))
  if (se && length(par_hat) > 0L && is.finite(value) && value < 1e11) {
    se_vec <- mp_try_se(obj, par_hat, nobs = length(x))
  }
  resid <- try(marpoly_residuals(x, causal = causal_hat, noncausal = noncausal_hat, demean = FALSE), silent = TRUE)
  resid_sd <- if (!inherits(resid, "try-error")) stats::sd(resid) else NA_real_
  df <- data.frame(
    allocation_id = allocation_row$allocation_id,
    model = allocation_row$model,
    r = r,
    s = s,
    objective = comps$value,
    convergence = convergence,
    residual_sd = resid_sd,
    stringsAsFactors = FALSE
  )
  mp_add_fit_lists(df, causal_hat, noncausal_hat, par_hat, start, se_vec, comps$kappas, message)
}

mp_add_fit_lists <- function(df, causal_hat, noncausal_hat, par_hat, start, se_vec, kappas, message) {
  df$causal <- list(as.numeric(Re(causal_hat)))
  df$noncausal <- list(as.numeric(Re(noncausal_hat)))
  df$parameters <- list(as.numeric(Re(par_hat)))
  df$start <- list(as.numeric(Re(start)))
  df$std_error <- list(as.numeric(se_vec))
  df$kappa <- list(kappas)
  df$optimizer_message <- list(message)
  df
}

mp_filter_allocations <- function(allocs, allocations) {
  if (identical(allocations, "all") || is.null(allocations)) {
    return(allocs)
  }
  if (is.character(allocations)) {
    keep <- allocs$model %in% allocations
    return(allocs[keep, , drop = FALSE])
  }
  if (is.matrix(allocations) || is.data.frame(allocations)) {
    rs <- as.data.frame(allocations)
    names(rs)[1:2] <- c("r", "s")
    keep <- rep(FALSE, nrow(allocs))
    for (i in seq_len(nrow(rs))) {
      keep <- keep | (allocs$r == rs$r[i] & allocs$s == rs$s[i])
    }
    return(allocs[keep, , drop = FALSE])
  }
  if (is.list(allocations)) {
    keep <- rep(FALSE, nrow(allocs))
    for (a in allocations) {
      keep <- keep | (allocs$r == a[1L] & allocs$s == a[2L])
    }
    return(allocs[keep, , drop = FALSE])
  }
  stop("allocations must be 'all', model labels, a two-column matrix/data frame, or a list of c(r, s) pairs.", call. = FALSE)
}

#' All nonempty combinations of spectral orders
#'
#' @param orders Candidate orders, usually `c(2, 3, 4)`.
#' @return A list of integer vectors.
#' @export
marpoly_order_sets <- function(orders = c(2, 3, 4)) {
  orders <- mp_parse_orders(orders)
  out <- list()
  for (m in seq_along(orders)) {
    cmb <- utils::combn(orders, m, simplify = FALSE)
    out <- c(out, cmb)
  }
  out
}

#' Estimate and rank MAR root allocations by polyspectral minimum distance
#'
#' This is the main estimation function. It first estimates a pseudo-causal AR(p)
#' representation, enumerates all feasible root allocations, estimates each
#' allocation with a frequency-domain objective based on selected polyspectral
#' orders, and ranks the allocations by the objective value.
#'
#' @param x Numeric time series.
#' @param p AR order. If NULL, selected by `marpoly_select_p`.
#' @param p_max Maximum AR order used when `p = NULL`.
#' @param orders Spectral orders used in the objective. Examples: `c(2,3)`,
#'   `c(2,4)`, or `c(2,3,4)`.
#' @param weights Optional objective weights, same length as `orders`.
#' @param grid_size Frequency grid size. Use a smaller value for order 4.
#' @param allocations Which root allocations to estimate. Default estimates all.
#' @param optimizer "single" for root-factorized starting values with local
#'   optimization, or "sann" for simulated annealing followed by a local polish.
#' @param maxit Maximum optimizer iterations.
#' @param criterion Information criterion for selecting p when `p = NULL`.
#' @param demean If TRUE, demean the series before estimation.
#' @param include_mean Include an intercept in the preliminary AR regression.
#' @param seed Optional seed for simulated annealing.
#' @param se If TRUE, compute numerical-Hessian standard errors.
#' @return Object of class `marpoly_fit`.
#' @export
marpoly_fit <- function(x, p = NULL, p_max = 6, orders = c(2, 3), weights = NULL,
                        grid_size = NULL, allocations = "all",
                        optimizer = c("single", "sann"), maxit = 3000,
                        criterion = c("BIC", "AIC"), demean = TRUE,
                        include_mean = FALSE, seed = NULL, se = TRUE) {
  optimizer <- match.arg(optimizer)
  criterion <- match.arg(criterion)
  orders <- mp_parse_orders(orders)
  if (is.null(weights)) {
    weights <- rep(1, length(orders))
  }
  weights <- as.numeric(weights)
  if (length(weights) != length(orders) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be a nonnegative numeric vector with the same length as orders.", call. = FALSE)
  }
  x_clean <- mp_clean_series(x, demean = demean)
  if (is.null(p)) {
    selection <- marpoly_select_p(x_clean, p_max = p_max, criterion = criterion, demean = FALSE, include_mean = include_mean)
    p <- selection$selected_p
    ar_fit <- selection$selected_fit
  } else {
    selection <- NULL
    ar_fit <- marpoly_ar_ols(x_clean, p = p, demean = FALSE, include_mean = include_mean)
  }
  eta <- ar_fit$ar
  allocs <- marpoly_root_allocations(eta)
  allocs <- mp_filter_allocations(allocs, allocations)
  if (nrow(allocs) == 0L) {
    stop("No root allocations selected.", call. = FALSE)
  }
  pdata <- mp_prepare_periodograms(x_clean, orders = orders, grid_size = grid_size, demean = FALSE)
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
  }
  rows <- vector("list", nrow(allocs))
  for (i in seq_len(nrow(allocs))) {
    rows[[i]] <- mp_fit_one_allocation(x_clean, allocs[i, ], pdata = pdata, eta = eta, orders = orders,
                                       weights = weights, optimizer = optimizer, maxit = maxit, se = se)
  }
  results <- do.call(rbind, rows)
  ord <- order(results$objective, na.last = TRUE)
  results <- results[ord, , drop = FALSE]
  results$rank <- seq_len(nrow(results))
  results <- results[, c("rank", setdiff(names(results), "rank"))]
  rownames(results) <- NULL
  out <- list(
    x = x_clean,
    p = p,
    orders = orders,
    weights = weights,
    optimizer = optimizer,
    preliminary_ar = ar_fit,
    order_selection = selection,
    allocations = allocs,
    results = results,
    selected_index = 1L,
    selected = results[1L, , drop = FALSE],
    periodograms = pdata,
    grid_size = grid_size,
    demean = demean
  )
  class(out) <- "marpoly_fit"
  out
}

#' Alias for the main MAR estimation function
#'
#' @inheritParams marpoly_fit
#' @export
marpoly_fit_all <- function(...) {
  marpoly_fit(...)
}

#' Rank table for a fitted MAR model
#'
#' @param fit A `marpoly_fit` object.
#' @param n Number of rows to return.
#' @return Data frame of ranked allocations.
#' @export
marpoly_rank <- function(fit, n = 10) {
  if (!inherits(fit, "marpoly_fit")) stop("fit must be a marpoly_fit object.", call. = FALSE)
  tab <- fit$results[, c("rank", "allocation_id", "model", "r", "s", "objective", "residual_sd", "convergence")]
  tab$causal <- vapply(fit$results$causal, mp_format_vector, character(1L))
  tab$noncausal <- vapply(fit$results$noncausal, mp_format_vector, character(1L))
  out <- utils::head(tab, n)
  class(out) <- c("marpoly_rank", class(out))
  out
}

#' @export
print.marpoly_rank <- function(x, ...) {
  NextMethod("print")
  invisible(x)
}

#' @export
coef.marpoly_fit <- function(object, allocation = "selected", ...) {
  if (identical(allocation, "selected")) {
    row <- object$results[object$selected_index, ]
  } else {
    id <- as.integer(allocation)
    pos <- which(object$results$allocation_id == id)
    if (length(pos) != 1L) stop("allocation id not found.", call. = FALSE)
    row <- object$results[pos, ]
  }
  out <- row$parameters[[1L]]
  names(out) <- c(if (row$r > 0L) paste0("causal", seq_len(row$r)) else character(0),
                  if (row$s > 0L) paste0("noncausal", seq_len(row$s)) else character(0))
  out
}

#' @export
summary.marpoly_fit <- function(object, competitors = 5, ...) {
  row <- object$results[object$selected_index, ]
  est <- row$parameters[[1L]]
  se <- row$std_error[[1L]]
  names <- c(if (row$r > 0L) paste0("causal", seq_len(row$r)) else character(0),
             if (row$s > 0L) paste0("noncausal", seq_len(row$s)) else character(0))
  if (length(est) == 0L) {
    coef_table <- data.frame()
  } else {
    stat <- est / se
    pval <- 2 * (1 - stats::pnorm(abs(stat)))
    coef_table <- data.frame(
      estimate = est,
      std.error = se,
      statistic = stat,
      p.value = pval,
      row.names = names
    )
  }
  out <- list(
    model = row$model,
    allocation_id = row$allocation_id,
    rank = row$rank,
    objective = row$objective,
    residual_sd = row$residual_sd,
    coefficients = coef_table,
    p = object$p,
    orders = object$orders,
    weights = object$weights,
    optimizer = object$optimizer,
    kappa = row$kappa[[1L]],
    competitors = marpoly_rank(object, n = competitors),
    preliminary_ar = object$preliminary_ar
  )
  class(out) <- "summary_marpoly_fit"
  out
}

#' @export
print.summary_marpoly_fit <- function(x, ...) {
  cat("marpoly frequency-domain MAR estimation\n")
  cat("  selected model: ", x$model, " (allocation ", x$allocation_id, ")\n", sep = "")
  cat("  pseudo-causal p: ", x$p, "\n", sep = "")
  cat("  spectral orders: ", paste(x$orders, collapse = ", "), "\n", sep = "")
  cat("  objective:       ", round(x$objective, 6), "\n", sep = "")
  cat("  residual sd:     ", round(x$residual_sd, 6), "\n", sep = "")
  if (nrow(x$coefficients) > 0L) {
    cat("\nCoefficients:\n")
    print(x$coefficients)
  }
  cat("\nEstimated cumulants used in the selected objective:\n")
  print(x$kappa)
  cat("\nRanked root allocations (main competitors):\n")
  print(x$competitors)
  invisible(x)
}

#' @export
print.marpoly_fit <- function(x, ...) {
  s <- summary(x)
  print(s)
  invisible(x)
}
