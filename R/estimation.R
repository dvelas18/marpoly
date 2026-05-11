mp_weights_to_R <- function(weights, orders) {
  if (is.null(weights)) return(NULL)
  weights <- as.numeric(weights)
  if (length(weights) == 1L) weights <- rep(weights, length(orders))
  if (length(weights) != length(orders)) {
    stop("weights must have length one or the same length as orders.", call. = FALSE)
  }
  out <- c(R2 = 0, R3 = 0, R4 = 0)
  out[paste0("R", orders)] <- weights
  out
}

mp_grid_control <- function(x, orders, grid_size, control) {
  # The core estimator uses domain-specific frequency rules from MAR_objective_prepare().
  # grid_size is kept only as a light compatibility shortcut for fourth-order slices;
  # it should not replace the bispectral/trispectral selection logic.
  if (is.null(grid_size)) return(control)
  n <- length(.mar_clean_series(x, demean = FALSE))
  g <- as.integer(grid_size[1L])
  if (!is.finite(g) || g < 1L) return(control)
  if (4L %in% orders && is.null(control$trispec_max_index)) {
    control$trispec_max_index <- min(g, floor((n - 1L) / 2L))
  }
  control
}

mp_specs_from_allocations <- function(allocations, p) {
  if (is.null(allocations) || identical(allocations, "all")) return(MAR_all_specs_from_p(p))
  all_specs <- MAR_all_specs_from_p(p)
  all_specs$model <- vapply(seq_len(nrow(all_specs)), function(i) .mar_model_label(all_specs$r[i], all_specs$s[i]), character(1))
  if (is.character(allocations)) {
    keep <- all_specs$model %in% allocations | paste0("MAR(", all_specs$r, ",", all_specs$s, ")") %in% allocations
    specs <- all_specs[keep, c("r", "s"), drop = FALSE]
  } else if (is.matrix(allocations) || is.data.frame(allocations)) {
    specs <- as.data.frame(allocations)
    names(specs)[1:2] <- c("r", "s")
    specs <- specs[, c("r", "s"), drop = FALSE]
  } else if (is.list(allocations)) {
    specs <- do.call(rbind, lapply(allocations, function(z) data.frame(r = as.integer(z[1L]), s = as.integer(z[2L]))))
  } else {
    stop("allocations must be 'all', model labels, a two-column matrix/data frame, or a list of c(r, s) pairs.", call. = FALSE)
  }
  if (nrow(specs) == 0L) stop("No MAR specifications matched allocations.", call. = FALSE)
  specs$r <- as.integer(specs$r)
  specs$s <- as.integer(specs$s)
  if (any(specs$r < 0L | specs$s < 0L | specs$r + specs$s != p)) {
    stop("Every selected allocation must satisfy r >= 0, s >= 0, and r + s = p.", call. = FALSE)
  }
  unique(specs)
}

mp_select_given_p <- function(y, p, eta_common, second_order, orders, specs, control) {
  fits <- vector("list", nrow(specs))
  root_tables <- vector("list", nrow(specs))
  base_seed <- control$random_seed %||% NULL
  for (i in seq_len(nrow(specs))) {
    r_i <- specs$r[i]
    s_i <- specs$s[i]
    control_i <- control
    control_i$eta_norm <- eta_common
    if (!is.null(base_seed)) control_i$random_seed <- base_seed + i - 1L
    fits[[i]] <- MAR_spectral_estimate(y = y, r = r_i, s = s_i, orders = orders, control = control_i)
    root_tables[[i]] <- MAR_root_combination_table(fits[[i]])
    root_tables[[i]]$spec_index <- i
  }
  names(fits) <- vapply(seq_len(nrow(specs)), function(i) .mar_model_label(specs$r[i], specs$s[i]), character(1))
  all_root_allocations <- .mar_bind_rows_fill(root_tables)
  if (nrow(all_root_allocations) > 0L) {
    all_root_allocations$selected_p <- p
    all_root_allocations <- .mar_reorder_cols(
      all_root_allocations,
      c("selected_p", "spec_index", "model", "r", "s", "orders", "root_combination",
        "start_source", "root_allocation", "value", "R2", "R3", "R4", "kappa2", "kappa3", "kappa4",
        "stationary", "stationarity_violation", "convergence", "n_freq", "n_pairs", "n_tuples")
    )
    all_root_allocations <- all_root_allocations[order(ifelse(is.finite(all_root_allocations$value),
                                                              all_root_allocations$value, Inf)), , drop = FALSE]
    rownames(all_root_allocations) <- NULL
  }
  report <- MAR_best_models_table(fits)
  selected <- report[1L, , drop = FALSE]
  out <- list(selected_p = p, eta_pseudo_causal = eta_common, second_order = second_order,
              specs = specs, all_root_allocations = all_root_allocations, report = report,
              best_by_model = report, selected = selected, fits = fits)
  class(out) <- c("marpoly_fit", "MAR_selection", "list")
  out
}

#' All nonempty combinations of spectral orders
#'
#' @param orders Candidate orders, usually `c(2, 3, 4)`.
#' @return A list of integer vectors.
#' @export
marpoly_order_sets <- function(orders = c(2, 3, 4)) {
  orders <- .mar_check_orders(orders)
  out <- list()
  for (m in seq_along(orders)) out <- c(out, utils::combn(orders, m, simplify = FALSE))
  out
}

#' Estimate and rank MAR root allocations by polyspectral minimum distance
#'
#' This wrapper uses the same frequency-domain objective and root-allocation
#' starting-value logic as `MAR_spectral_estimate()` and
#' `MAR_select_all_by_second_order()`.
#'
#' @param x Numeric time series.
#' @param p AR order. If NULL, selected by second-order BIC/AIC/AICc.
#' @param p_max Maximum AR order used when `p = NULL`.
#' @param orders Spectral orders used in the objective.
#' @param weights Optional objective weights in the same order as `orders`.
#' @param grid_size Convenience control for coarser frequency grids.
#' @param allocations Which MAR(r,s) specifications to estimate.
#' @param optimizer "single" for root-factorized starts and Nelder-Mead, or
#'   "sann" for simulated annealing.
#' @param maxit Maximum optimizer iterations.
#' @param criterion Information criterion for selecting p when `p = NULL`.
#' @param demean If TRUE, demean the input series.
#' @param include_mean Kept for backward compatibility; the spectral workflow
#'   demeans through `demean`.
#' @param seed Optional seed.
#' @param se If TRUE, compute sandwich standard errors.
#' @param control Additional control passed to `MAR_spectral_estimate()`.
#' @param ar_method Method passed to `stats::ar()` for preliminary AR estimates.
#' @param ... Additional control arguments.
#' @return Object of class `marpoly_fit`.
#' @export
marpoly_fit <- function(x, p = NULL, p_max = 6, orders = c(2, 3), weights = NULL,
                        grid_size = NULL, allocations = "all",
                        optimizer = c("single", "sann"), maxit = 3000,
                        criterion = c("BIC", "AIC", "AICc"), demean = TRUE,
                        include_mean = FALSE, seed = NULL, se = TRUE,
                        control = list(), ar_method = c("yule-walker", "burg", "ols", "mle"), ...) {
  optimizer <- match.arg(optimizer)
  criterion <- match.arg(criterion)
  ar_method <- match.arg(ar_method)
  orders <- .mar_check_orders(orders)
  extra <- list(...)
  if (length(extra) > 0L) control <- utils::modifyList(control, extra)
  control <- mp_grid_control(x, orders, grid_size, control)
  control$weights <- mp_weights_to_R(weights, orders) %||% control$weights
  control$demean <- demean
  control$eta_method <- ar_method
  control$compute_se <- se
  if (is.null(control$compute_errors)) control$compute_errors <- TRUE
  if (is.null(control$random_seed)) control$random_seed <- seed
  if (optimizer == "sann") {
    control$optim_method <- "SANN"
    control$optim_control <- utils::modifyList(list(maxit = maxit), control$optim_control %||% list())
  } else {
    control$optim_method <- "Nelder-Mead"
    control$optim_control <- utils::modifyList(list(maxit = maxit, reltol = 1e-8), control$optim_control %||% list())
  }
  y_clean <- .mar_clean_series(x, demean = demean)
  if (is.null(p)) {
    second_order <- MAR_select_p_second_order(y_clean, max_p = p_max, ic = tolower(criterion),
                                              method = ar_method, demean = FALSE)
    p <- second_order$p
    eta <- second_order$eta
  } else {
    p <- .mar_check_nonnegative_integer(p, "p")
    eta <- MAR_estimate_causal_eta(y_clean, p = p, method = ar_method, demean = FALSE)
    second_order <- list(p = p, eta = eta, criterion = "fixed", method = ar_method,
                         table = data.frame(p = p, sigma2 = NA_real_, ic = NA_real_, converged = TRUE),
                         ljung_box = data.frame())
  }
  specs <- mp_specs_from_allocations(allocations, p)
  out <- mp_select_given_p(y = y_clean, p = p, eta_common = eta, second_order = second_order,
                           orders = orders, specs = specs, control = control)
  out$x <- y_clean
  out$orders <- orders
  out$weights <- control$weights
  out$optimizer <- optimizer
  out$grid_size <- grid_size
  out$demean <- demean
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
#' @param root_allocations If TRUE, rank every local root allocation; otherwise
#'   rank the best result by MAR(r,s) model.
#' @return Data frame of ranked models or allocations.
#' @export
marpoly_rank <- function(fit, n = 10, root_allocations = FALSE) {
  if (!inherits(fit, "marpoly_fit") && !inherits(fit, "MAR_selection")) {
    stop("fit must be a marpoly_fit or MAR_selection object.", call. = FALSE)
  }
  tab <- if (isTRUE(root_allocations)) fit$all_root_allocations else fit$report
  if (is.null(tab) || nrow(tab) == 0L) return(data.frame())
  tab <- tab[order(ifelse(is.finite(tab$value), tab$value, Inf)), , drop = FALSE]
  tab$rank <- seq_len(nrow(tab))
  tab$objective <- tab$value
  tab <- .mar_reorder_cols(tab, c("rank", "model", "r", "s", "orders", "objective", "R2", "R3", "R4", "kappa2", "kappa3", "kappa4", "convergence", "parameters"))
  out <- utils::head(tab, n)
  class(out) <- c("marpoly_rank", class(out))
  out
}

#' @export
print.marpoly_rank <- function(x, ...) {
  NextMethod("print")
  invisible(x)
}

mp_selected_fit <- function(object) {
  if (!inherits(object, "MAR_selection") && !inherits(object, "marpoly_fit")) {
    stop("object must be a MAR selection object.", call. = FALSE)
  }
  nm <- object$selected$model[1]
  fit <- object$fits[[nm]]
  if (is.null(fit)) fit <- object$fits[[1L]]
  fit
}

#' @export
coef.marpoly_fit <- function(object, allocation = "selected", ...) {
  if (identical(allocation, "selected")) return(mp_selected_fit(object)$theta)
  id <- as.integer(allocation)
  tab <- object$all_root_allocations
  pos <- which(tab$root_combination == id)
  if (length(pos) < 1L) stop("allocation id not found.", call. = FALSE)
  row <- tab[pos[1L], , drop = FALSE]
  names <- .mar_make_names(row$r[1L], row$s[1L])
  out <- as.numeric(row[1L, names])
  names(out) <- names
  out
}

#' @export
summary.marpoly_fit <- function(object, competitors = 5, ...) {
  fit <- mp_selected_fit(object)
  out <- list(
    model = object$selected$model[1],
    selected_p = object$selected_p,
    orders = fit$orders,
    objective = fit$value,
    R2 = fit$R2,
    R3 = fit$R3,
    R4 = fit$R4,
    kappa2 = fit$kappa2,
    kappa3 = fit$kappa3,
    kappa4 = fit$kappa4,
    coefficients = fit$parameter_table,
    residual_sd = if (!is.null(fit$errors)) stats::sd(fit$errors) else NA_real_,
    convergence = fit$convergence,
    se_method = fit$se_method,
    competitors = marpoly_rank(object, n = competitors),
    root_allocations = marpoly_rank(object, n = competitors, root_allocations = TRUE),
    second_order = object$second_order
  )
  class(out) <- "summary_marpoly_fit"
  out
}

#' @export
print.summary_marpoly_fit <- function(x, ...) {
  cat("marpoly frequency-domain MAR estimation\n")
  cat("  selected model: ", x$model, "\n", sep = "")
  cat("  pseudo-causal p: ", x$selected_p, "\n", sep = "")
  cat("  spectral orders: ", paste(x$orders, collapse = ", "), "\n", sep = "")
  cat("  objective:       ", round(x$objective, 6), "\n", sep = "")
  cat("  R2/R3/R4:        ", paste(round(c(x$R2, x$R3, x$R4), 6), collapse = " / "), "\n", sep = "")
  cat("  residual sd:     ", round(x$residual_sd, 6), "\n", sep = "")
  if (!is.null(x$coefficients) && nrow(x$coefficients) > 0L) {
    cat("\nCoefficients:\n")
    print(x$coefficients, row.names = FALSE)
  }
  cat("\nMain competing MAR(r,s) models:\n")
  print(x$competitors)
  cat("\nMain root-allocation/local-solution competitors:\n")
  print(x$root_allocations)
  invisible(x)
}

#' @export
print.marpoly_fit <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' Estimate all MAR(r,s) specifications for a fixed pseudo-causal order p
#'
#' @param y Numeric time series.
#' @param p Fixed pseudo-causal AR order.
#' @param orders Spectral orders used in the objective.
#' @param ar_method Preliminary causal AR estimation method.
#' @param specs Optional data frame with columns r and s.
#' @param control Control list passed to the spectral estimator.
#' @param allocations Optional specification filter.
#' @param ... Additional control arguments.
#' @return A MAR selection object.
#' @export
MAR_select_all_fixed_p <- function(y, p, orders = c(2, 3),
                                   ar_method = c("yule-walker", "burg", "ols", "mle"),
                                   specs = NULL, control = list(), allocations = "all", ...) {
  ar_method <- match.arg(ar_method)
  orders <- .mar_check_orders(orders)
  p <- .mar_check_nonnegative_integer(p, "p")
  extra <- list(...)
  if (length(extra) > 0L) control <- utils::modifyList(control, extra)
  y_clean <- .mar_clean_series(y, demean = control$demean %||% TRUE)
  eta <- MAR_estimate_causal_eta(y_clean, p = p, method = ar_method, demean = FALSE)
  second_order <- list(p = p, eta = eta, criterion = "fixed", method = ar_method,
                       table = data.frame(p = p, sigma2 = NA_real_, ic = NA_real_, converged = TRUE),
                       ljung_box = data.frame())
  if (is.null(specs)) specs <- mp_specs_from_allocations(allocations, p)
  out <- mp_select_given_p(y = y_clean, p = p, eta_common = eta, second_order = second_order,
                           orders = orders, specs = specs, control = control)
  out$x <- y_clean
  class(out) <- c("marpoly_fit", "MAR_selection", "list")
  out
}

#' @export
print.MAR_selection <- function(x, ...) {
  if (!inherits(x, "marpoly_fit")) class(x) <- c("marpoly_fit", class(x))
  print.marpoly_fit(x, ...)
}

#' @export
print.MAR_fit <- function(x, ...) {
  cat("MAR spectral estimate\n")
  cat("  model: ", .mar_model_label(x$r, x$s), "\n", sep = "")
  cat("  spectral orders: ", paste(x$orders, collapse = ", "), "\n", sep = "")
  cat("  objective: ", round(x$value, 6), "\n", sep = "")
  if (!is.null(x$parameter_table) && nrow(x$parameter_table) > 0L) {
    cat("\nCoefficients:\n")
    print(x$parameter_table, row.names = FALSE)
  }
  invisible(x)
}
