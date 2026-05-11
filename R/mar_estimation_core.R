# ============================================================
# Spectral / bispectral / trispectral minimum-distance
# estimation for MAR(r,s)
#
# Orders:
#   orders = c(2)       spectrum only
#   orders = c(3)       bispectrum only
#   orders = c(4)       trispectrum only
#   orders = c(2,3)     spectrum + bispectrum
#   orders = c(2,4)     spectrum + trispectrum
#   orders = c(3,4)     bispectrum + trispectrum
#   orders = c(2,3,4)   all components
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.mar_check_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 0 || x != floor(x)) {
    stop(paste0(name, " must be a non-negative integer."), call. = FALSE)
  }
  as.integer(x)
}

.mar_check_positive_integer <- function(x, name, min_value = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < min_value || x != floor(x)) {
    stop(paste0(name, " must be an integer >= ", min_value, "."), call. = FALSE)
  }
  as.integer(x)
}

.mar_check_orders <- function(orders = c(2, 3), order = NULL) {
  if (!is.null(order)) orders <- order
  if (is.character(orders)) {
    if (length(orders) == 1L && tolower(orders) == "all") return(c(2L, 3L, 4L))
    stop("orders must be numeric, for example c(2,3,4), or 'all'.", call. = FALSE)
  }
  if (!is.numeric(orders) || length(orders) == 0L ||
      any(!is.finite(orders)) || any(orders != floor(orders))) {
    stop("orders must contain integers from 2 to 4.", call. = FALSE)
  }
  orders <- sort(unique(as.integer(orders)))
  if (!all(orders %in% c(2L, 3L, 4L))) {
    stop("orders can only contain 2, 3, and/or 4.", call. = FALSE)
  }
  orders
}

.mar_clean_series <- function(y, demean = TRUE) {
  if (!is.numeric(y) && !is.complex(y)) stop("y must be numeric or complex.", call. = FALSE)
  y <- as.vector(y)
  if (length(y) < 5L) stop("y must contain at least five observations.", call. = FALSE)
  if (any(!is.finite(Re(y))) || any(!is.finite(Im(y)))) {
    stop("y contains non-finite values.", call. = FALSE)
  }
  if (demean) y <- y - mean(y)
  y
}

.mar_clean_coefs <- function(x, name, zero_tol = 1e-12) {
  if (is.null(x)) return(numeric(0))
  if (!is.numeric(x) && !is.complex(x)) {
    stop(paste0(name, " must be numeric, complex, or NULL."), call. = FALSE)
  }
  x <- as.vector(x)
  if (length(x) == 0L) return(numeric(0))
  if (any(!is.finite(Re(x))) || any(!is.finite(Im(x)))) {
    stop(paste0(name, " contains non-finite values."), call. = FALSE)
  }
  while (length(x) > 0L && Mod(x[length(x)]) <= zero_tol) x <- x[-length(x)]
  if (length(x) > 0L && max(abs(Im(x))) <= 1e-10) x <- Re(x)
  x
}

.mar_split_theta <- function(theta, r, s) {
  theta <- as.numeric(theta)
  if (length(theta) != r + s) {
    stop(paste0("theta must have length r + s = ", r + s, "."), call. = FALSE)
  }
  list(
    phi_plus = if (r > 0L) theta[seq_len(r)] else numeric(0),
    phi_star = if (s > 0L) theta[r + seq_len(s)] else numeric(0)
  )
}

.mar_make_names <- function(r, s) {
  c(
    if (r > 0L) paste0("phi_plus_", seq_len(r)) else character(0),
    if (s > 0L) paste0("phi_star_", seq_len(s)) else character(0)
  )
}

.mar_model_label <- function(r, s) {
  if (s == 0L) return(paste0("AR(", r, ",0) causal"))
  if (r == 0L) return(paste0("AR(0,", s, ") noncausal"))
  paste0("MAR(", r, ",", s, ")")
}

.mar_poly_multiply <- function(a, b) {
  out <- rep(0 + 0i, length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    out[i:(i + length(b) - 1L)] <- out[i:(i + length(b) - 1L)] + a[i] * b
  }
  out
}

.mar_roots <- function(coefs, zero_tol = 1e-12) {
  coefs <- as.vector(coefs)
  if (length(coefs) == 0L) return(complex(0))
  p <- c(1, -coefs)
  while (length(p) > 1L && Mod(p[length(p)]) <= zero_tol) p <- p[-length(p)]
  if (length(p) <= 1L) return(complex(0))
  base::polyroot(p)
}

.mar_stationarity_info <- function(phi_plus, phi_star, root_tol = 1e-8, zero_tol = 1e-12) {
  roots_plus <- .mar_roots(phi_plus, zero_tol = zero_tol)
  roots_star <- .mar_roots(phi_star, zero_tol = zero_tol)
  viol_plus <- if (length(roots_plus) == 0L) 0 else sum(pmax(0, 1 + root_tol - Mod(roots_plus))^2)
  viol_star <- if (length(roots_star) == 0L) 0 else sum(pmax(0, 1 + root_tol - Mod(roots_star))^2)
  list(
    stationary = (viol_plus + viol_star) == 0,
    violation = viol_plus + viol_star,
    roots_phi_plus = roots_plus,
    roots_phi_star = roots_star,
    min_root_modulus_phi_plus = if (length(roots_plus) == 0L) Inf else min(Mod(roots_plus)),
    min_root_modulus_phi_star = if (length(roots_star) == 0L) Inf else min(Mod(roots_star))
  )
}

.mar_dft <- function(y) {
  T <- length(y)
  omega_all <- 2 * pi * (0:(T - 1L)) / T
  stats::fft(y) * exp(-1i * omega_all)
}

.mar_transfer <- function(omega, phi_plus = NULL, phi_star = NULL) {
  phi_plus <- .mar_clean_coefs(phi_plus, "phi_plus")
  phi_star <- .mar_clean_coefs(phi_star, "phi_star")
  dims <- dim(omega)
  omega_vec <- as.vector(omega)
  if (!is.numeric(omega_vec) || any(!is.finite(omega_vec))) stop("omega must be numeric and finite.", call. = FALSE)
  z <- exp(-1i * omega_vec)
  denom <- rep(1 + 0i, length(z))
  if (length(phi_plus) > 0L) {
    powers_plus <- outer(z, seq_along(phi_plus), function(a, b) a^b)
    denom <- denom * (1 - as.vector(powers_plus %*% phi_plus))
  }
  if (length(phi_star) > 0L) {
    powers_star <- outer(z, -seq_along(phi_star), function(a, b) a^b)
    denom <- denom * (1 - as.vector(powers_star %*% phi_star))
  }
  psi <- 1 / denom
  if (!is.null(dims)) psi <- array(psi, dim = dims)
  psi
}

.mar_bind_rows_fill <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(data.frame())
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows2 <- lapply(rows, function(x) {
    missing_names <- setdiff(all_names, names(x))
    for (nm in missing_names) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, rows2)
}

.mar_reorder_cols <- function(x, first_cols) {
  first_cols <- intersect(first_cols, names(x))
  x[, c(first_cols, setdiff(names(x), first_cols)), drop = FALSE]
}

.mar_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(default)
  x <- as.numeric(x[1L])
  if (!is.finite(x)) return(default)
  x
}

.mar_pinv <- function(A, tol = NULL) {
  s <- svd(A)
  if (is.null(tol)) tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
  keep <- s$d > tol
  if (!any(keep)) return(matrix(0, nrow = ncol(A), ncol = nrow(A)))
  s$v[, keep, drop = FALSE] %*% diag(1 / s$d[keep], nrow = sum(keep)) %*% t(s$u[, keep, drop = FALSE])
}

# ============================================================
# 1. Frequency domains
# ============================================================

MAR_spectrum_indices <- function(T, step = 1L) {
  T <- .mar_check_positive_integer(T, "T", min_value = 5L)
  step <- .mar_check_positive_integer(step, "step", min_value = 1L)
  seq.int(1L, floor((T - 1L) / 2L), by = step)
}

MAR_bispectrum_pairs <- function(T, step = NULL, max_pairs = 2e5, include_boundary = TRUE) {
  T <- .mar_check_positive_integer(T, "T", min_value = 5L)
  if (is.null(step)) {
    approx_pairs <- T^2 / 24
    step <- max(1L, as.integer(ceiling(sqrt(approx_pairs / max_pairs))))
  }
  step <- .mar_check_positive_integer(step, "step", min_value = 1L)
  k1_max <- floor((T - 1L) / 2L)
  out1 <- integer(0)
  out2 <- integer(0)
  for (k1 in seq.int(1L, k1_max, by = step)) {
    if (include_boundary) k2_max <- min(k1, T - 2L * k1) else k2_max <- min(k1 - 1L, T - 2L * k1 - 1L)
    if (k2_max >= 1L) {
      k2_grid <- seq.int(1L, k2_max, by = step)
      out1 <- c(out1, rep(k1, length(k2_grid)))
      out2 <- c(out2, k2_grid)
    }
  }
  cbind(k1 = out1, k2 = out2)
}

.center_index <- function(k, T) {
  mod <- as.integer(k %% T)
  out <- mod
  out[mod > T / 2] <- mod[mod > T / 2] - T
  as.integer(out)
}

.is_nyquist <- function(k, T) {
  if (T %% 2L != 0L) return(rep(FALSE, length(k)))
  as.integer(k %% T) == T / 2L
}

.lex_less <- function(a, b) {
  for (i in seq_along(a)) {
    if (a[i] < b[i]) return(TRUE)
    if (a[i] > b[i]) return(FALSE)
  }
  FALSE
}

.canonical_quad <- function(q, T) {
  q1 <- sort(.center_index(q, T))
  q2 <- sort(.center_index(-q, T))
  if (.lex_less(q2, q1)) q2 else q1
}

MAR_trispectrum_tuples <- function(T, step = NULL, max_tuples = 50000,
                                   max_index = NULL, candidate_index = NULL,
                                   include_negative = TRUE, exclude_pairwise = TRUE,
                                   unique = TRUE) {
  T <- .mar_check_positive_integer(T, "T", min_value = 8L)
  max_tuples <- .mar_check_positive_integer(max_tuples, "max_tuples", min_value = 1L)
  if (!is.null(candidate_index)) {
    candidates <- as.integer(candidate_index)
    candidates <- unique(.center_index(candidates, T))
    candidates <- candidates[candidates != 0L]
    candidates <- candidates[!.is_nyquist(candidates, T)]
    if (include_negative) {
      candidates <- unique(c(candidates, -candidates))
      candidates <- candidates[candidates != 0L]
      candidates <- candidates[!.is_nyquist(candidates, T)]
    }
    candidates <- sort(unique(candidates))
  } else {
    if (is.null(max_index)) max_index <- min(40L, floor((T - 1L) / 2L))
    max_index <- .mar_check_positive_integer(max_index, "max_index", min_value = 1L)
    max_index <- min(max_index, floor((T - 1L) / 2L))
    if (is.null(step)) {
      target_candidates <- max(8L, floor((2 * max_tuples)^(1 / 3)))
      step <- max(1L, as.integer(ceiling((2 * max_index) / target_candidates)))
    }
    step <- .mar_check_positive_integer(step, "step", min_value = 1L)
    positive <- seq.int(1L, max_index, by = step)
    candidates <- if (include_negative) sort(unique(c(-rev(positive), positive))) else positive
    candidates <- candidates[!.is_nyquist(candidates, T)]
  }
  if (length(candidates) == 0L) stop("No trispectrum candidate frequencies were selected.", call. = FALSE)
  rows <- matrix(NA_integer_, nrow = max_tuples, ncol = 3L)
  colnames(rows) <- c("k1", "k2", "k3")
  seen <- new.env(parent = emptyenv())
  count <- 0L
  hit_cap <- FALSE
  for (a in candidates) {
    for (b in candidates) {
      for (c in candidates) {
        d <- -(a + b + c)
        q <- .center_index(c(a, b, c, d), T)
        if (any(q == 0L)) next
        if (any(.is_nyquist(q, T))) next
        if (exclude_pairwise) {
          pair_sums <- c(q[1L] + q[2L], q[1L] + q[3L], q[1L] + q[4L],
                         q[2L] + q[3L], q[2L] + q[4L], q[3L] + q[4L])
          if (any((pair_sums %% T) == 0L)) next
        }
        q_can <- if (unique) .canonical_quad(q, T) else q
        key <- paste(q_can, collapse = ",")
        if (unique && exists(key, envir = seen, inherits = FALSE)) next
        if (unique) assign(key, TRUE, envir = seen)
        count <- count + 1L
        rows[count, ] <- q_can[1:3]
        if (count >= max_tuples) {
          hit_cap <- TRUE
          break
        }
      }
      if (count >= max_tuples) break
    }
    if (count >= max_tuples) break
  }
  if (count == 0L) {
    stop("No trispectrum tuples remained after applying uniqueness and degeneracy filters.", call. = FALSE)
  }
  out <- rows[seq_len(count), , drop = FALSE]
  attr(out, "hit_cap") <- hit_cap
  attr(out, "rule") <- paste0(
    "Canonical fourth-order frequency tuples. Quadruple is (k1,k2,k3,-k1-k2-k3). ",
    "Permutations and conjugate duplicates are removed. Pairwise-zero frequency pairs are excluded."
  )
  if (hit_cap) {
    warning(paste0("The trispectrum tuple cap was reached: ", max_tuples,
                   ". Increase max_tuples or use a coarser trispec_step/trispec_max_index if needed."),
            call. = FALSE)
  }
  out
}

# ============================================================
# 2. Second-order pseudo-causal AR(p) selection
# ============================================================

.mar_ar_residuals <- function(y, eta) {
  eta <- as.numeric(eta)
  p <- length(eta)
  T <- length(y)
  if (p == 0L) return(as.numeric(y))
  e <- rep(NA_real_, T)
  for (t in (p + 1L):T) e[t] <- y[t] - sum(eta * y[t - seq_len(p)])
  e
}

MAR_estimate_causal_eta <- function(y, p, method = c("yule-walker", "burg", "ols", "mle"), demean = TRUE) {
  p <- .mar_check_nonnegative_integer(p, "p")
  if (p == 0L) return(numeric(0))
  y <- .mar_clean_series(y, demean = demean)
  method <- match.arg(method)
  fit <- stats::ar(y, order.max = p, aic = FALSE, method = method, demean = FALSE)
  eta <- as.numeric(fit$ar)
  if (length(eta) != p || any(!is.finite(eta))) {
    stop("Could not estimate the preliminary causal AR coefficients eta.", call. = FALSE)
  }
  eta
}

MAR_select_p_second_order <- function(y, max_p = 8L, ic = c("bic", "aic", "aicc"),
                                      method = c("yule-walker", "burg", "ols", "mle"),
                                      demean = TRUE) {
  max_p <- .mar_check_nonnegative_integer(max_p, "max_p")
  ic <- match.arg(tolower(ic[1L]), c("bic", "aic", "aicc"))
  method <- match.arg(method)
  y <- .mar_clean_series(y, demean = demean)
  T <- length(y)
  rows <- vector("list", max_p + 1L)
  fits <- vector("list", max_p + 1L)
  for (p in 0:max_p) {
    if (p == 0L) {
      eta <- numeric(0)
      residuals <- y
      sigma2 <- mean(residuals^2)
      n_eff <- T
      converged <- TRUE
    } else {
      fit_p <- tryCatch(stats::ar(y, order.max = p, aic = FALSE, method = method, demean = FALSE),
                        error = function(e) NULL)
      if (is.null(fit_p) || length(fit_p$ar) != p || any(!is.finite(fit_p$ar))) {
        eta <- rep(NA_real_, p)
        residuals <- rep(NA_real_, T)
        sigma2 <- Inf
        n_eff <- T - p
        converged <- FALSE
      } else {
        eta <- as.numeric(fit_p$ar)
        residuals <- .mar_ar_residuals(y, eta)
        residuals_used <- residuals[is.finite(residuals)]
        sigma2 <- mean(residuals_used^2)
        n_eff <- length(residuals_used)
        converged <- TRUE
      }
    }
    penalty <- switch(ic, aic = 2 * p / n_eff, bic = log(n_eff) * p / n_eff,
                      aicc = 2 * p / max(1L, n_eff - p - 1L))
    ic_value <- log(sigma2) + penalty
    rows[[p + 1L]] <- data.frame(p = p, sigma2 = sigma2, ic = ic_value, converged = converged)
    fits[[p + 1L]] <- list(p = p, eta = eta, residuals = residuals, sigma2 = sigma2,
                           ic = ic_value, converged = converged)
  }
  table <- do.call(rbind, rows)
  best_idx <- which.min(ifelse(is.finite(table$ic), table$ic, Inf))
  best <- fits[[best_idx]]
  residuals <- best$residuals
  residuals <- residuals[is.finite(residuals)]
  T_res <- length(residuals)
  lb_long_lag <- min(max(10L, best$p + 1L), max(1L, floor(T_res / 4L)))
  lb_lags <- unique(c(1L, 2L, lb_long_lag))
  ljung_box <- data.frame(
    lag = lb_lags,
    p_value = vapply(lb_lags, function(lag) {
      if (lag <= best$p || length(residuals) <= lag + best$p) return(NA_real_)
      tryCatch(stats::Box.test(residuals, lag = lag, type = "Ljung-Box", fitdf = best$p)$p.value,
               error = function(e) NA_real_)
    }, numeric(1))
  )
  list(p = best$p, eta = best$eta, residuals = best$residuals, sigma2 = best$sigma2,
       criterion = ic, method = method, table = table, ljung_box = ljung_box)
}

MAR_all_specs_from_p <- function(p) {
  p <- .mar_check_nonnegative_integer(p, "p")
  data.frame(r = p:0, s = 0:p)
}

# ============================================================
# 3. Root-allocation starts
# ============================================================

.mar_coefs_from_roots <- function(roots, imag_tol = 1e-7) {
  roots <- as.vector(roots)
  if (length(roots) == 0L) return(numeric(0))
  poly <- c(1 + 0i)
  for (root in roots) poly <- .mar_poly_multiply(poly, c(1, -1 / root))
  coefs <- -poly[-1L]
  if (length(coefs) > 0L && max(abs(Im(coefs))) <= imag_tol) coefs <- Re(coefs)
  coefs
}

.mar_coefs_from_lambdas <- function(lambda, imag_tol = 1e-7) {
  lambda <- as.vector(lambda)
  if (length(lambda) == 0L) return(numeric(0))
  poly <- c(1 + 0i)
  for (l in lambda) poly <- .mar_poly_multiply(poly, c(1, -l))
  coefs <- -poly[-1L]
  if (length(coefs) > 0L && max(abs(Im(coefs))) <= imag_tol) coefs <- Re(coefs)
  coefs
}

.mar_root_units <- function(roots, imag_tol = 1e-7) {
  n <- length(roots)
  if (n == 0L) return(list())
  unused <- rep(TRUE, n)
  units <- list()
  for (i in seq_len(n)) {
    if (!unused[i]) next
    unused[i] <- FALSE
    if (abs(Im(roots[i])) <= imag_tol) {
      units[[length(units) + 1L]] <- i
    } else {
      candidates <- which(unused)
      if (length(candidates) == 0L) {
        units[[length(units) + 1L]] <- i
      } else {
        distances <- Mod(roots[candidates] - Conj(roots[i]))
        j <- candidates[which.min(distances)]
        if (Mod(roots[j] - Conj(roots[i])) <= 100 * imag_tol) {
          unused[j] <- FALSE
          units[[length(units) + 1L]] <- c(i, j)
        } else {
          units[[length(units) + 1L]] <- i
        }
      }
    }
  }
  units
}

.mar_unit_combinations <- function(unit_sizes, target) {
  out <- list()
  rec <- function(i, chosen, total) {
    if (total == target) {
      out[[length(out) + 1L]] <<- chosen
      return()
    }
    if (total > target || i > length(unit_sizes)) return()
    rec(i + 1L, chosen, total)
    rec(i + 1L, c(chosen, i), total + unit_sizes[i])
  }
  rec(1L, integer(0), 0L)
  out
}

MAR_auto_starts_from_eta <- function(eta, r, s, max_starts = Inf, imag_tol = 1e-7,
                                     zero_tol = 1e-12, return_allocations = FALSE) {
  r <- .mar_check_nonnegative_integer(r, "r")
  s <- .mar_check_nonnegative_integer(s, "s")
  p <- r + s
  eta <- as.numeric(eta)
  if (length(eta) != p) stop("eta must have length r + s.", call. = FALSE)
  if (p == 0L) {
    if (return_allocations) return(list(list(theta = numeric(0), source = "automatic", allocation = "MAR(0,0): no roots")))
    return(list(numeric(0)))
  }
  roots <- .mar_roots(eta, zero_tol = zero_tol)
  units <- .mar_root_units(roots, imag_tol = imag_tol)
  unit_sizes <- lengths(units)
  combos <- .mar_unit_combinations(unit_sizes, target = r)
  if (length(combos) == 0L) {
    warning(paste0("No real-valued root allocation was found for MAR(", r, ",", s,"). Using a zero start instead."),
            call. = FALSE)
    if (return_allocations) return(list(list(theta = rep(0, p), source = "automatic-zero", allocation = "No admissible real root allocation")))
    return(list(rep(0, p)))
  }
  if (is.finite(max_starts) && length(combos) > max_starts) combos <- combos[seq_len(max_starts)]
  starts <- list()
  all_indices <- seq_along(roots)
  for (combo in combos) {
    plus_idx <- if (length(combo) == 0L) integer(0) else unlist(units[combo], use.names = FALSE)
    star_idx <- setdiff(all_indices, plus_idx)
    phi_plus <- .mar_coefs_from_roots(roots[plus_idx], imag_tol = imag_tol)
    phi_star <- .mar_coefs_from_roots(roots[star_idx], imag_tol = imag_tol)
    theta <- c(phi_plus, phi_star)
    if (length(theta) == p && max(abs(Im(theta))) <= imag_tol) {
      theta <- as.numeric(Re(theta))
      allocation <- paste0("causal root indices = {", paste(plus_idx, collapse = ","),
                           "}; noncausal root indices = {", paste(star_idx, collapse = ","), "}")
      starts[[length(starts) + 1L]] <- if (return_allocations) {
        list(theta = theta, source = "automatic-root-allocation", allocation = allocation)
      } else {
        theta
      }
    }
  }
  if (length(starts) == 0L) {
    warning("All automatic starts were complex. Using a zero start instead.", call. = FALSE)
    if (return_allocations) return(list(list(theta = rep(0, p), source = "automatic-zero", allocation = "All automatic starts were complex")))
    return(list(rep(0, p)))
  }
  theta_mat <- do.call(rbind, lapply(starts, function(x) if (return_allocations) x$theta else x))
  keep <- !duplicated(round(theta_mat, 12))
  starts[keep]
}

.mar_random_starts <- function(r, s, n, radius = 0.8, seed = NULL) {
  if (n <= 0L) return(list())
  if (!is.null(seed)) set.seed(seed)
  starts <- vector("list", n)
  for (i in seq_len(n)) {
    lambda_plus <- if (r > 0L) stats::runif(r, -radius, radius) else numeric(0)
    lambda_star <- if (s > 0L) stats::runif(s, -radius, radius) else numeric(0)
    theta <- c(.mar_coefs_from_lambdas(lambda_plus), .mar_coefs_from_lambdas(lambda_star))
    starts[[i]] <- list(theta = as.numeric(theta), source = "random", allocation = "random stationary coefficient start")
  }
  starts
}

.mar_normalize_starts <- function(start, p) {
  if (is.null(start)) return(list())
  if (is.vector(start) && !is.list(start)) {
    start <- as.numeric(start)
    if (length(start) != p) stop("start must have length r + s.", call. = FALSE)
    return(list(list(theta = start, source = "manual", allocation = "manual start")))
  }
  if (is.list(start) && !(is.data.frame(start))) {
    out <- lapply(start, as.numeric)
    bad <- which(vapply(out, length, integer(1)) != p)
    if (length(bad) > 0L) stop("Every start vector must have length r + s.", call. = FALSE)
    return(lapply(seq_along(out), function(i) list(theta = out[[i]], source = "manual", allocation = paste0("manual start ", i))))
  }
  if (is.matrix(start) || is.data.frame(start)) {
    start <- as.matrix(start)
    if (ncol(start) != p) stop("start matrix/data.frame must have r + s columns.", call. = FALSE)
    return(lapply(seq_len(nrow(start)), function(i) list(theta = as.numeric(start[i, ]), source = "manual", allocation = paste0("manual row ", i))))
  }
  stop("start must be a vector, list, matrix, data.frame, or NULL.", call. = FALSE)
}

# ============================================================
# 4. Estimated model errors / innovations
# ============================================================

.mar_laurent_operator <- function(phi_plus, phi_star) {
  phi_plus <- .mar_clean_coefs(phi_plus, "phi_plus")
  phi_star <- .mar_clean_coefs(phi_star, "phi_star")
  r <- length(phi_plus)
  s <- length(phi_star)
  a_powers <- c(0L, if (r > 0L) seq_len(r) else integer(0))
  a_coef <- c(1, if (r > 0L) -phi_plus else numeric(0))
  b_powers <- c(0L, if (s > 0L) -seq_len(s) else integer(0))
  b_coef <- c(1, if (s > 0L) -phi_star else numeric(0))
  powers <- seq.int(-s, r)
  out <- setNames(rep(0, length(powers)), as.character(powers))
  for (i in seq_along(a_powers)) {
    for (j in seq_along(b_powers)) {
      k <- a_powers[i] + b_powers[j]
      out[as.character(k)] <- out[as.character(k)] + a_coef[i] * b_coef[j]
    }
  }
  out
}

MAR_model_errors <- function(y, phi_plus = NULL, phi_star = NULL, demean = TRUE) {
  y <- .mar_clean_series(y, demean = demean)
  phi_plus <- .mar_clean_coefs(phi_plus, "phi_plus")
  phi_star <- .mar_clean_coefs(phi_star, "phi_star")
  r <- length(phi_plus)
  s <- length(phi_star)
  T <- length(y)
  if (T <= r + s + 1L) stop("The series is too short to compute MAR errors.", call. = FALSE)
  op <- .mar_laurent_operator(phi_plus, phi_star)
  powers <- as.integer(names(op))
  t_index <- seq.int(r + 1L, T - s)
  eps_hat <- vapply(t_index, function(t) sum(op * y[t - powers]), numeric(1))
  names(eps_hat) <- as.character(t_index)
  eps_hat
}

MAR_error_diagnostics <- function(errors, lags = c(1L, 2L, 5L, 10L), powers = 1:4) {
  e <- as.numeric(errors)
  e <- e[is.finite(e)]
  if (length(e) < 10L) return(list(summary = data.frame(), ljung_box = data.frame()))
  m <- mean(e)
  sd_e <- stats::sd(e)
  summary <- data.frame(n = length(e), mean = m, sd = sd_e,
                        skewness = if (sd_e > 0) mean((e - m)^3) / sd_e^3 else NA_real_,
                        kurtosis = if (sd_e > 0) mean((e - m)^4) / sd_e^4 else NA_real_)
  rows <- list()
  for (pow in powers) {
    x <- e^pow
    x <- x - mean(x)
    for (lag in lags) {
      pval <- if (length(x) > lag + 2L) {
        tryCatch(stats::Box.test(x, lag = lag, type = "Ljung-Box", fitdf = 0)$p.value,
                 error = function(err) NA_real_)
      } else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(power = pow, lag = lag, p_value = pval)
    }
  }
  list(summary = summary, ljung_box = do.call(rbind, rows))
}

# ============================================================
# 5. Objective preparation
# ============================================================

.mar_parse_weights <- function(weights, orders) {
  if (is.null(weights)) weights <- c(R2 = 1, R3 = 1, R4 = 1)
  if (length(weights) == 1L) weights <- c(R2 = weights[1], R3 = weights[1], R4 = weights[1])
  if (is.null(names(weights))) {
    tmp <- c(R2 = 0, R3 = 0, R4 = 0)
    tmp[seq_along(weights)] <- weights
    weights <- tmp
  } else {
    tmp <- c(R2 = 0, R3 = 0, R4 = 0)
    valid <- intersect(names(weights), names(tmp))
    tmp[valid] <- weights[valid]
    if (!all(names(weights) %in% names(tmp))) {
      unnamed <- weights[!(names(weights) %in% names(tmp))]
      if (length(unnamed) > 0L) tmp[seq_along(unnamed)] <- unnamed
    }
    weights <- tmp
  }
  if (any(!is.finite(weights)) || any(weights < 0)) stop("weights must be finite and non-negative.", call. = FALSE)
  if (!(2L %in% orders)) weights["R2"] <- 0
  if (!(3L %in% orders)) weights["R3"] <- 0
  if (!(4L %in% orders)) weights["R4"] <- 0
  weights
}

MAR_objective_prepare <- function(y, r, s, orders = c(2, 3), order = NULL,
                                  freq_index = NULL, pairs = NULL, tuples = NULL,
                                  freq_step = 1L, bispec_step = NULL,
                                  max_bispec_pairs = 2e5, include_bispec_boundary = TRUE,
                                  trispec_step = NULL, max_trispec_tuples = 50000,
                                  trispec_max_index = NULL, trispec_candidate_index = NULL,
                                  demean = TRUE, eta_norm = NULL,
                                  eta_method = c("yule-walker", "burg", "ols", "mle"),
                                  weights = NULL, scaling = c("mean", "paper"),
                                  check_stationarity = TRUE, root_tol = 1e-8,
                                  zero_tol = 1e-12, denom_tol = 1e-10,
                                  penalty = 1e12) {
  r <- .mar_check_nonnegative_integer(r, "r")
  s <- .mar_check_nonnegative_integer(s, "s")
  p <- r + s
  orders <- .mar_check_orders(orders = orders, order = order)
  scaling <- match.arg(scaling)
  eta_method <- match.arg(eta_method)
  weights <- .mar_parse_weights(weights, orders)
  y <- .mar_clean_series(y, demean = demean)
  T <- length(y)

  if (2L %in% orders) {
    if (is.null(freq_index)) idx2 <- MAR_spectrum_indices(T, step = freq_step) else {
      idx2 <- as.integer(freq_index)
      idx2 <- idx2[idx2 %% T != 0L]
      if (length(idx2) == 0L) stop("freq_index is empty after excluding zero frequency.", call. = FALSE)
    }
  } else idx2 <- integer(0)

  if (3L %in% orders) {
    if (is.null(pairs)) {
      pairs <- MAR_bispectrum_pairs(T = T, step = bispec_step, max_pairs = max_bispec_pairs,
                                    include_boundary = include_bispec_boundary)
    } else {
      pairs <- as.matrix(pairs)
      if (ncol(pairs) != 2L) stop("pairs must have two columns: k1 and k2.", call. = FALSE)
      pairs <- matrix(as.integer(pairs), ncol = 2L)
      colnames(pairs) <- c("k1", "k2")
    }
    if (nrow(pairs) == 0L) stop("No bispectrum frequency pairs were selected.", call. = FALSE)
    k1 <- as.integer(pairs[, 1L])
    k2 <- as.integer(pairs[, 2L])
    k3 <- -(k1 + k2)
  } else {
    pairs <- NULL
    k1 <- k2 <- k3 <- integer(0)
  }

  if (4L %in% orders) {
    if (is.null(tuples)) {
      tuples <- MAR_trispectrum_tuples(T = T, step = trispec_step, max_tuples = max_trispec_tuples,
                                       max_index = trispec_max_index,
                                       candidate_index = trispec_candidate_index,
                                       include_negative = TRUE, exclude_pairwise = TRUE, unique = TRUE)
    } else {
      tuples <- as.matrix(tuples)
      if (ncol(tuples) != 3L) stop("tuples must have three columns: k1, k2, and k3.", call. = FALSE)
      tuples <- matrix(as.integer(tuples), ncol = 3L)
      colnames(tuples) <- c("k1", "k2", "k3")
    }
    if (nrow(tuples) == 0L) stop("No trispectrum frequency tuples were selected.", call. = FALSE)
    h1 <- as.integer(tuples[, 1L]); h2 <- as.integer(tuples[, 2L]); h3 <- as.integer(tuples[, 3L]); h4 <- -(h1 + h2 + h3)
  } else {
    tuples <- NULL
    h1 <- h2 <- h3 <- h4 <- integer(0)
  }

  d_all <- .mar_dft(y)
  d_at <- function(idx) d_all[as.integer(idx %% T) + 1L]
  I2 <- if (2L %in% orders) d_at(idx2) * d_at(-idx2) / (2 * pi * T) else complex(0)
  I3 <- if (3L %in% orders) d_at(k1) * d_at(k2) * d_at(k3) / ((2 * pi)^2 * T) else complex(0)
  I4 <- if (4L %in% orders) d_at(h1) * d_at(h2) * d_at(h3) * d_at(h4) / ((2 * pi)^3 * T) else complex(0)

  if (is.null(eta_norm)) {
    eta_norm <- MAR_estimate_causal_eta(y = y, p = p, method = eta_method, demean = FALSE)
  } else {
    eta_norm <- as.numeric(eta_norm)
  }
  if (length(eta_norm) != p) stop("eta_norm must have length r + s.", call. = FALSE)

  idx_needed <- unique(c(idx2, -idx2, k1, k2, k3, h1, h2, h3, h4))
  idx_needed <- idx_needed[is.finite(idx_needed)]
  if (length(idx_needed) == 0L) stop("No frequency indices were selected.", call. = FALSE)
  idx_needed_mod <- sort(unique(as.integer(idx_needed %% T)))
  omega_needed <- 2 * pi * idx_needed_mod / T
  z_needed <- exp(-1i * omega_needed)
  n_freq <- length(idx_needed_mod)
  make_Z <- function(power_seq) {
    if (length(power_seq) == 0L) return(matrix(0 + 0i, nrow = n_freq, ncol = 0L))
    outer(z_needed, power_seq, function(a, b) a^b)
  }
  Z_plus <- make_Z(if (r > 0L) seq_len(r) else integer(0))
  Z_star <- make_Z(if (s > 0L) -seq_len(s) else integer(0))
  Z_eta <- make_Z(if (p > 0L) seq_len(p) else integer(0))
  idx_map <- rep(NA_integer_, T)
  idx_map[idx_needed_mod + 1L] <- seq_along(idx_needed_mod)
  pos <- function(idx) {
    out <- idx_map[as.integer(idx %% T) + 1L]
    if (any(is.na(out))) stop("Internal frequency-index mapping error.", call. = FALSE)
    out
  }
  eval_psi_theta <- function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) != p || any(!is.finite(theta))) return(NULL)
    parts <- .mar_split_theta(theta, r, s)
    denom <- rep(1 + 0i, n_freq)
    if (r > 0L) denom <- denom * (1 - as.vector(Z_plus %*% parts$phi_plus))
    if (s > 0L) denom <- denom * (1 - as.vector(Z_star %*% parts$phi_star))
    if (any(!is.finite(Re(denom))) || any(!is.finite(Im(denom))) || any(Mod(denom) < denom_tol)) return(NULL)
    1 / denom
  }
  eval_psi_eta <- function() {
    denom <- rep(1 + 0i, n_freq)
    if (p > 0L) denom <- denom * (1 - as.vector(Z_eta %*% eta_norm))
    if (any(Mod(denom) < denom_tol)) stop("The preliminary causal normalization has a near-zero denominator.", call. = FALSE)
    1 / denom
  }
  psi_eta <- eval_psi_eta()

  if (2L %in% orders) {
    pos2 <- pos(idx2); pos2_neg <- pos(-idx2)
    eta_prod2 <- psi_eta[pos2] * psi_eta[pos2_neg]
    den2 <- Mod(eta_prod2)^2
    if (any(den2 < denom_tol)) stop("Second-order normalization denominator is too small.", call. = FALSE)
  } else { pos2 <- pos2_neg <- integer(0); den2 <- numeric(0) }

  if (3L %in% orders) {
    pos31 <- pos(k1); pos32 <- pos(k2); pos33 <- pos(k3)
    eta_prod3 <- psi_eta[pos31] * psi_eta[pos32] * psi_eta[pos33]
    den3 <- Mod(eta_prod3)^2
    if (any(den3 < denom_tol)) stop("Third-order normalization denominator is too small.", call. = FALSE)
  } else { pos31 <- pos32 <- pos33 <- integer(0); den3 <- numeric(0) }

  if (4L %in% orders) {
    pos41 <- pos(h1); pos42 <- pos(h2); pos43 <- pos(h3); pos44 <- pos(h4)
    eta_prod4 <- psi_eta[pos41] * psi_eta[pos42] * psi_eta[pos43] * psi_eta[pos44]
    den4 <- Mod(eta_prod4)^2
    if (any(den4 < denom_tol)) stop("Fourth-order normalization denominator is too small.", call. = FALSE)
  } else { pos41 <- pos42 <- pos43 <- pos44 <- integer(0); den4 <- numeric(0) }

  evaluate <- function(theta, details = TRUE) {
    theta <- as.numeric(theta)
    if (length(theta) != p || any(!is.finite(theta))) return(if (details) list(value = penalty) else penalty)
    parts <- .mar_split_theta(theta, r, s)
    stat <- .mar_stationarity_info(parts$phi_plus, parts$phi_star, root_tol = root_tol, zero_tol = zero_tol)
    if (check_stationarity && !stat$stationary) {
      val <- penalty * (1 + stat$violation)
      return(if (details) list(value = val, R2 = NA_real_, R3 = NA_real_, R4 = NA_real_,
                               kappa2 = NA_real_, kappa3 = NA_real_, kappa4 = NA_real_,
                               stationary = FALSE, stationarity_violation = stat$violation,
                               phi_plus = parts$phi_plus, phi_star = parts$phi_star) else val)
    }
    psi <- eval_psi_theta(theta)
    if (is.null(psi)) return(if (details) list(value = penalty) else penalty)
    R2 <- R3 <- R4 <- NA_real_
    kappa2 <- kappa3 <- kappa4 <- NA_real_

    if (2L %in% orders) {
      psi_prod2 <- psi[pos2] * psi[pos2_neg]
      if (any(Mod(psi_prod2) < denom_tol)) return(if (details) list(value = penalty) else penalty)
      kappa2 <- 2 * pi * mean(Re(I2 / psi_prod2))
      if (!is.finite(kappa2) || kappa2 <= 0) return(if (details) list(value = penalty) else penalty)
      S2 <- kappa2 / (2 * pi) * psi_prod2
      raw2 <- Mod(I2 - S2)^2 / den2
      R2 <- if (scaling == "mean") mean(raw2) else sum(raw2) / (4 * T)
    }
    if (3L %in% orders) {
      psi_prod3 <- psi[pos31] * psi[pos32] * psi[pos33]
      if (any(Mod(psi_prod3) < denom_tol)) return(if (details) list(value = penalty) else penalty)
      kappa3 <- (2 * pi)^2 * mean(Re(I3 / psi_prod3))
      if (!is.finite(kappa3)) return(if (details) list(value = penalty) else penalty)
      S3 <- kappa3 / ((2 * pi)^2) * psi_prod3
      raw3 <- Mod(I3 - S3)^2 / den3
      R3 <- if (scaling == "mean") mean(raw3) else (2 * pi) * sum(raw3) / (6 * T^2)
    }
    if (4L %in% orders) {
      psi_prod4 <- psi[pos41] * psi[pos42] * psi[pos43] * psi[pos44]
      if (any(Mod(psi_prod4) < denom_tol)) return(if (details) list(value = penalty) else penalty)
      kappa4 <- (2 * pi)^3 * mean(Re(I4 / psi_prod4))
      if (!is.finite(kappa4)) return(if (details) list(value = penalty) else penalty)
      S4 <- kappa4 / ((2 * pi)^3) * psi_prod4
      raw4 <- Mod(I4 - S4)^2 / den4
      R4 <- if (scaling == "mean") mean(raw4) else (2 * pi)^2 * sum(raw4) / (8 * T^3)
    }
    value <- weights["R2"] * ifelse(is.na(R2), 0, R2) +
      weights["R3"] * ifelse(is.na(R3), 0, R3) +
      weights["R4"] * ifelse(is.na(R4), 0, R4)
    if (!is.finite(value)) value <- penalty
    if (!details) return(value)
    list(value = value, R2 = R2, R3 = R3, R4 = R4,
         kappa2 = kappa2, kappa3 = kappa3, kappa4 = kappa4,
         stationary = stat$stationary, stationarity_violation = stat$violation,
         phi_plus = parts$phi_plus, phi_star = parts$phi_star)
  }

  moment_vector <- function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) != p || any(!is.finite(theta))) return(NULL)
    parts <- .mar_split_theta(theta, r, s)
    stat <- .mar_stationarity_info(parts$phi_plus, parts$phi_star, root_tol = root_tol, zero_tol = zero_tol)
    if (check_stationarity && !stat$stationary) return(NULL)
    psi <- eval_psi_theta(theta)
    if (is.null(psi)) return(NULL)
    out <- numeric(0)
    if (2L %in% orders) {
      psi_prod2 <- psi[pos2] * psi[pos2_neg]
      if (any(Mod(psi_prod2) < denom_tol)) return(NULL)
      kappa2 <- 2 * pi * mean(Re(I2 / psi_prod2))
      if (!is.finite(kappa2) || kappa2 <= 0) return(NULL)
      S2 <- kappa2 / (2 * pi) * psi_prod2
      e2 <- (I2 - S2) / sqrt(den2)
      out <- c(out, sqrt(weights["R2"]) * Re(e2))
    }
    if (3L %in% orders) {
      psi_prod3 <- psi[pos31] * psi[pos32] * psi[pos33]
      if (any(Mod(psi_prod3) < denom_tol)) return(NULL)
      kappa3 <- (2 * pi)^2 * mean(Re(I3 / psi_prod3))
      if (!is.finite(kappa3)) return(NULL)
      S3 <- kappa3 / ((2 * pi)^2) * psi_prod3
      e3 <- (I3 - S3) / sqrt(den3)
      out <- c(out, sqrt(weights["R3"]) * Re(e3), sqrt(weights["R3"]) * Im(e3))
    }
    if (4L %in% orders) {
      psi_prod4 <- psi[pos41] * psi[pos42] * psi[pos43] * psi[pos44]
      if (any(Mod(psi_prod4) < denom_tol)) return(NULL)
      kappa4 <- (2 * pi)^3 * mean(Re(I4 / psi_prod4))
      if (!is.finite(kappa4)) return(NULL)
      S4 <- kappa4 / ((2 * pi)^3) * psi_prod4
      e4 <- (I4 - S4) / sqrt(den4)
      out <- c(out, sqrt(weights["R4"]) * Re(e4), sqrt(weights["R4"]) * Im(e4))
    }
    out[is.finite(out)]
  }

  objective <- function(theta) evaluate(theta, details = FALSE)
  list(objective = objective, evaluate = evaluate, moment_vector = moment_vector,
       y = y, T = T, r = r, s = s, p = p, orders = orders,
       freq_index = idx2, pairs = pairs, tuples = tuples, eta_norm = eta_norm,
       I2 = I2, I3 = I3, I4 = I4, weights = weights, scaling = scaling,
       n_freq = length(idx2), n_pairs = if (is.null(pairs)) 0L else nrow(pairs),
       n_tuples = if (is.null(tuples)) 0L else nrow(tuples),
       convention = paste0("R2 uses periodogram/spectrum; R3 uses biperiodogram/bispectrum; ",
                           "R4 uses triperiodogram/trispectral density on canonical non-degenerate tuples."))
}

# ============================================================
# 6. Standard errors
# ============================================================

MAR_standard_errors <- function(setup, theta, fd_eps = 1e-5) {
  theta <- as.numeric(theta)
  p <- length(theta)
  if (p == 0L) return(list(se = numeric(0), vcov = matrix(0, 0, 0), method = "empty parameter vector"))
  q0 <- setup$moment_vector(theta)
  if (is.null(q0) || length(q0) <= p) {
    return(list(se = rep(NA_real_, p), vcov = matrix(NA_real_, p, p),
                method = "sandwich failed: invalid moment vector"))
  }
  n_mom <- length(q0)
  D <- matrix(NA_real_, nrow = n_mom, ncol = p)
  for (j in seq_len(p)) {
    h <- fd_eps * max(1, abs(theta[j]))
    th_plus <- theta; th_minus <- theta
    th_plus[j] <- th_plus[j] + h
    th_minus[j] <- th_minus[j] - h
    q_plus <- setup$moment_vector(th_plus)
    q_minus <- setup$moment_vector(th_minus)
    if (is.null(q_plus) || is.null(q_minus) || length(q_plus) != n_mom || length(q_minus) != n_mom) {
      return(list(se = rep(NA_real_, p), vcov = matrix(NA_real_, p, p),
                  method = "sandwich failed: invalid finite-difference step"))
    }
    D[, j] <- (q_plus - q_minus) / (2 * h)
  }
  A <- crossprod(D) / n_mom
  B <- crossprod(D, D * (q0^2)) / n_mom
  A_inv <- .mar_pinv(A)
  vcov <- A_inv %*% B %*% A_inv / setup$T
  vcov <- (vcov + t(vcov)) / 2
  se <- sqrt(pmax(0, diag(vcov)))
  list(se = se, vcov = vcov, method = "finite-difference sandwich based on selected spectral moments")
}

# ============================================================
# 7. Main estimation function
# ============================================================

MAR_spectral_estimate <- function(y, r, s, orders = c(2, 3), order = NULL,
                                  start = NULL, control = list(), ...) {
  r <- .mar_check_nonnegative_integer(r, "r")
  s <- .mar_check_nonnegative_integer(s, "s")
  p <- r + s
  orders <- .mar_check_orders(orders = orders, order = order)
  extra <- list(...)
  default_control <- list(
    auto_start = TRUE, max_auto_starts = Inf, n_random_starts = 0L,
    random_start_radius = 0.8, random_seed = NULL,
    optim_method = "Nelder-Mead", optim_control = list(maxit = 3000, reltol = 1e-8),
    compute_errors = FALSE, compute_objective_errors = FALSE,
    compute_se = TRUE, return_setup = FALSE
  )
  control <- utils::modifyList(default_control, control)
  if (length(extra) > 0L) control <- utils::modifyList(control, extra)
  prepare_names <- c("freq_index", "pairs", "tuples", "freq_step", "bispec_step", "max_bispec_pairs",
                     "include_bispec_boundary", "trispec_step", "max_trispec_tuples", "trispec_max_index",
                     "trispec_candidate_index", "demean", "eta_norm", "eta_method", "weights", "scaling",
                     "check_stationarity", "root_tol", "zero_tol", "denom_tol", "penalty")
  prepare_args <- control[intersect(names(control), prepare_names)]
  setup <- do.call(MAR_objective_prepare, c(list(y = y, r = r, s = s, orders = orders), prepare_args))
  names_par <- .mar_make_names(r, s)

  if (p == 0L) {
    details <- setup$evaluate(numeric(0), details = TRUE)
    out <- list(theta = numeric(0), phi_plus = numeric(0), phi_star = numeric(0),
                r = r, s = s, p = p, orders = orders, value = details$value,
                R2 = details$R2, R3 = details$R3, R4 = details$R4,
                kappa2 = details$kappa2, kappa3 = details$kappa3, kappa4 = details$kappa4,
                stationary = TRUE, stationarity_violation = 0, eta_norm = setup$eta_norm,
                convergence = 0, message = "No parameters to estimate: MAR(0,0).",
                parameter_table = data.frame(), all_results = list())
    if (isTRUE(control$compute_errors)) {
      out$errors <- MAR_model_errors(y = setup$y, phi_plus = numeric(0), phi_star = numeric(0), demean = FALSE)
      out$error_diagnostics <- MAR_error_diagnostics(out$errors)
    }
    if (isTRUE(control$compute_objective_errors)) out$objective_errors <- setup$moment_vector(numeric(0))
    if (isTRUE(control$return_setup)) out$setup <- setup
    class(out) <- c("MAR_fit", "list")
    return(out)
  }

  starts <- list()
  if (isTRUE(control$auto_start)) {
    starts <- c(starts, MAR_auto_starts_from_eta(eta = setup$eta_norm, r = r, s = s,
                                                 max_starts = control$max_auto_starts,
                                                 return_allocations = TRUE))
  }
  starts <- c(starts, .mar_normalize_starts(start, p))
  n_random_starts <- as.integer(control$n_random_starts)
  if (n_random_starts > 0L) {
    starts <- c(starts, .mar_random_starts(r = r, s = s, n = n_random_starts,
                                           radius = control$random_start_radius,
                                           seed = control$random_seed))
  }
  if (length(starts) == 0L) starts <- list(list(theta = rep(0, p), source = "zero", allocation = "zero start"))
  start_mat <- do.call(rbind, lapply(starts, function(x) x$theta))
  keep <- !duplicated(round(start_mat, 12))
  starts <- starts[keep]

  results <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    st <- as.numeric(starts[[i]]$theta)
    opt <- tryCatch(stats::optim(par = st, fn = setup$objective, method = control$optim_method,
                                 control = control$optim_control),
                    error = function(e) list(par = st, value = Inf, convergence = 999,
                                             message = conditionMessage(e)))
    details <- setup$evaluate(opt$par, details = TRUE)
    results[[i]] <- list(start = st, source = starts[[i]]$source, allocation = starts[[i]]$allocation,
                         optim = opt, details = details, value = details$value)
  }
  values <- vapply(results, function(x) x$value, numeric(1))
  best_idx <- which.min(ifelse(is.finite(values), values, Inf))
  best <- results[[best_idx]]
  theta_hat <- as.numeric(best$optim$par)
  names(theta_hat) <- names_par
  parts <- .mar_split_theta(theta_hat, r, s)
  se_info <- if (isTRUE(control$compute_se)) MAR_standard_errors(setup, theta_hat) else {
    list(se = rep(NA_real_, length(theta_hat)), vcov = matrix(NA_real_, length(theta_hat), length(theta_hat)), method = "not computed")
  }
  names(se_info$se) <- names_par
  rownames(se_info$vcov) <- names_par
  colnames(se_info$vcov) <- names_par
  t_value <- as.numeric(theta_hat) / as.numeric(se_info$se)
  parameter_table <- data.frame(parameter = names_par, estimate = as.numeric(theta_hat),
                                std_error = as.numeric(se_info$se), t_value = t_value,
                                p_value = 2 * (1 - stats::pnorm(abs(t_value))))

  out <- list(theta = theta_hat, phi_plus = parts$phi_plus, phi_star = parts$phi_star,
              r = r, s = s, p = p, orders = orders, value = best$details$value,
              R2 = best$details$R2, R3 = best$details$R3, R4 = best$details$R4,
              kappa2 = best$details$kappa2, kappa3 = best$details$kappa3, kappa4 = best$details$kappa4,
              stationary = best$details$stationary,
              stationarity_violation = best$details$stationarity_violation,
              roots_phi_plus = .mar_roots(parts$phi_plus), roots_phi_star = .mar_roots(parts$phi_star),
              eta_norm = setup$eta_norm, n_freq = setup$n_freq, n_pairs = setup$n_pairs,
              n_tuples = setup$n_tuples, best_start_index = best_idx, best_start_source = best$source,
              best_root_allocation = best$allocation, convergence = best$optim$convergence,
              message = if ("message" %in% names(best$optim)) best$optim$message else NULL,
              parameter_table = parameter_table, std_errors = se_info$se,
              vcov = se_info$vcov, se_method = se_info$method, all_results = results)
  if (isTRUE(control$compute_errors)) {
    out$errors <- MAR_model_errors(y = setup$y, phi_plus = parts$phi_plus, phi_star = parts$phi_star, demean = FALSE)
    out$error_diagnostics <- MAR_error_diagnostics(out$errors)
  }
  if (isTRUE(control$compute_objective_errors)) out$objective_errors <- setup$moment_vector(theta_hat)
  if (isTRUE(control$return_setup)) out$setup <- setup
  class(out) <- c("MAR_fit", "list")
  out
}

# ============================================================
# 8. Tables and clean reports
# ============================================================

MAR_root_combination_table <- function(fit) {
  r <- fit$r; s <- fit$s; p <- r + s
  model_label <- .mar_model_label(r, s)
  par_names <- .mar_make_names(r, s)
  if (p == 0L || is.null(fit$all_results) || length(fit$all_results) == 0L) {
    out <- data.frame(model = model_label, r = r, s = s, orders = paste(fit$orders, collapse = "+"),
                      root_combination = 1L, start_source = "none", root_allocation = "none",
                      value = fit$value, R2 = fit$R2, R3 = fit$R3, R4 = fit$R4,
                      kappa2 = fit$kappa2, kappa3 = fit$kappa3, kappa4 = fit$kappa4,
                      stationary = fit$stationary %||% TRUE,
                      stationarity_violation = fit$stationarity_violation %||% 0,
                      convergence = fit$convergence, n_freq = fit$n_freq %||% NA_integer_,
                      n_pairs = fit$n_pairs %||% NA_integer_, n_tuples = fit$n_tuples %||% NA_integer_)
    return(out)
  }
  rows <- vector("list", length(fit$all_results))
  for (i in seq_along(fit$all_results)) {
    res <- fit$all_results[[i]]
    opt <- res$optim
    details <- res$details %||% list()
    theta_hat <- as.numeric(opt$par)
    start <- as.numeric(res$start)
    if (length(theta_hat) != p) theta_hat <- rep(NA_real_, p)
    if (length(start) != p) start <- rep(NA_real_, p)
    names(theta_hat) <- par_names
    names(start) <- par_names
    parts <- tryCatch(.mar_split_theta(theta_hat, r, s),
                      error = function(e) list(phi_plus = rep(NA_real_, r), phi_star = rep(NA_real_, s)))
    stat <- tryCatch(.mar_stationarity_info(parts$phi_plus, parts$phi_star),
                     error = function(e) list(stationary = NA, violation = NA_real_,
                                              min_root_modulus_phi_plus = NA_real_,
                                              min_root_modulus_phi_star = NA_real_))
    out <- data.frame(model = model_label, r = r, s = s, orders = paste(fit$orders, collapse = "+"),
                      root_combination = i, start_source = res$source %||% NA_character_,
                      root_allocation = res$allocation %||% NA_character_,
                      value = .mar_scalar(details$value, res$value), R2 = .mar_scalar(details$R2),
                      R3 = .mar_scalar(details$R3), R4 = .mar_scalar(details$R4),
                      kappa2 = .mar_scalar(details$kappa2), kappa3 = .mar_scalar(details$kappa3),
                      kappa4 = .mar_scalar(details$kappa4), stationary = stat$stationary,
                      stationarity_violation = .mar_scalar(stat$violation),
                      convergence = .mar_scalar(opt$convergence),
                      min_root_modulus_phi_plus = .mar_scalar(stat$min_root_modulus_phi_plus, Inf),
                      min_root_modulus_phi_star = .mar_scalar(stat$min_root_modulus_phi_star, Inf),
                      n_freq = fit$n_freq %||% NA_integer_, n_pairs = fit$n_pairs %||% NA_integer_,
                      n_tuples = fit$n_tuples %||% NA_integer_)
    for (j in seq_along(theta_hat)) {
      out[[par_names[j]]] <- theta_hat[j]
      out[[paste0("start_", par_names[j])]] <- start[j]
    }
    rows[[i]] <- out
  }
  .mar_bind_rows_fill(rows)
}

.mar_parameter_string <- function(parameter_table) {
  if (is.null(parameter_table) || nrow(parameter_table) == 0L) return("(none)")
  parts <- character(nrow(parameter_table))
  for (i in seq_len(nrow(parameter_table))) {
    est <- parameter_table$estimate[i]
    se <- parameter_table$std_error[i]
    if (is.finite(se)) {
      parts[i] <- paste0(parameter_table$parameter[i], " = ", sprintf("%.5f", est),
                         " (se ", sprintf("%.5f", se), ")")
    } else {
      parts[i] <- paste0(parameter_table$parameter[i], " = ", sprintf("%.5f", est))
    }
  }
  paste(parts, collapse = "; ")
}

MAR_best_models_table <- function(fits) {
  rows <- vector("list", length(fits))
  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    rows[[i]] <- data.frame(model = .mar_model_label(fit$r, fit$s), r = fit$r, s = fit$s,
                            p = fit$p, orders = paste(fit$orders, collapse = "+"),
                            value = fit$value, R2 = fit$R2, R3 = fit$R3, R4 = fit$R4,
                            kappa2 = fit$kappa2, kappa3 = fit$kappa3, kappa4 = fit$kappa4,
                            n_freq = fit$n_freq, n_pairs = fit$n_pairs, n_tuples = fit$n_tuples,
                            convergence = fit$convergence,
                            best_start_source = fit$best_start_source %||% NA_character_,
                            best_root_allocation = fit$best_root_allocation %||% NA_character_,
                            parameters = .mar_parameter_string(fit$parameter_table))
  }
  out <- .mar_bind_rows_fill(rows)
  out <- out[order(ifelse(is.finite(out$value), out$value, Inf)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

MAR_report <- function(selection, top = Inf) {
  if (!"report" %in% names(selection)) stop("selection does not contain a report table.", call. = FALSE)
  out <- selection$report
  if (is.finite(top)) out <- head(out, top)
  out
}

# ============================================================
# 9. Full logical selection: second-order p -> MAR(r,s)
# ============================================================

MAR_select_all_by_second_order <- function(y, orders = c(2, 3), order = NULL,
                                           max_p = 8L, ic = c("bic", "aic", "aicc"),
                                           ar_method = c("yule-walker", "burg", "ols", "mle"),
                                           specs = NULL, control = list(), ...) {
  orders <- .mar_check_orders(orders = orders, order = order)
  ic <- match.arg(tolower(ic[1L]), c("bic", "aic", "aicc"))
  ar_method <- match.arg(ar_method)
  extra <- list(...)
  if (length(extra) > 0L) control <- utils::modifyList(control, extra)
  second_order <- MAR_select_p_second_order(y = y, max_p = max_p, ic = ic,
                                            method = ar_method, demean = control$demean %||% TRUE)
  p <- second_order$p
  eta_common <- second_order$eta
  if (is.null(specs)) {
    specs <- MAR_all_specs_from_p(p)
  } else {
    specs <- as.data.frame(specs)
    if (!all(c("r", "s") %in% names(specs))) stop("specs must contain columns named r and s.", call. = FALSE)
    specs$r <- as.integer(specs$r); specs$s <- as.integer(specs$s)
    if (any(specs$r + specs$s != p)) stop("Every specification must satisfy r + s = selected p.", call. = FALSE)
  }
  fits <- vector("list", nrow(specs))
  root_tables <- vector("list", nrow(specs))
  base_seed <- control$random_seed %||% NULL
  for (i in seq_len(nrow(specs))) {
    r_i <- specs$r[i]; s_i <- specs$s[i]
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
    all_root_allocations <- .mar_reorder_cols(all_root_allocations,
      c("selected_p", "spec_index", "model", "r", "s", "orders", "root_combination",
        "start_source", "root_allocation", "value", "R2", "R3", "R4", "kappa2", "kappa3", "kappa4",
        "stationary", "stationarity_violation", "convergence", "n_freq", "n_pairs", "n_tuples"))
    all_root_allocations <- all_root_allocations[order(ifelse(is.finite(all_root_allocations$value),
                                                              all_root_allocations$value, Inf)), , drop = FALSE]
    rownames(all_root_allocations) <- NULL
  }
  report <- MAR_best_models_table(fits)
  selected <- report[1L, , drop = FALSE]
  out <- list(x = .mar_clean_series(y, demean = control$demean %||% TRUE),
              selected_p = p, eta_pseudo_causal = eta_common, second_order = second_order,
              specs = specs, all_root_allocations = all_root_allocations, report = report,
              best_by_model = report, selected = selected, fits = fits)
  class(out) <- c("MAR_selection", "list")
  out
}

MAR_select_all_by_auto_arima <- function(...) {
  MAR_select_all_by_second_order(...)
}

# ============================================================
# 10. Compare all possible order combinations
# ============================================================

MAR_all_order_sets <- function() {
  list(c(2), c(3), c(4), c(2, 3), c(2, 4), c(3, 4), c(2, 3, 4))
}

MAR_select_by_order_sets <- function(y, order_sets = MAR_all_order_sets(), max_p = 8L,
                                     ic = "bic", ar_method = "yule-walker", control = list(), ...) {
  extra <- list(...)
  if (length(extra) > 0L) control <- utils::modifyList(control, extra)
  selections <- vector("list", length(order_sets))
  rows <- vector("list", length(order_sets))
  for (i in seq_along(order_sets)) {
    ord_i <- .mar_check_orders(order_sets[[i]])
    selections[[i]] <- MAR_select_all_by_second_order(y = y, orders = ord_i, max_p = max_p,
                                                      ic = ic, ar_method = ar_method, control = control)
    sel <- selections[[i]]$selected
    rows[[i]] <- data.frame(order_set = paste(ord_i, collapse = "+"),
                            selected_p = selections[[i]]$selected_p,
                            selected_model = sel$model[1], r = sel$r[1], s = sel$s[1],
                            value = sel$value[1], R2 = sel$R2[1], R3 = sel$R3[1], R4 = sel$R4[1],
                            kappa2 = sel$kappa2[1], kappa3 = sel$kappa3[1], kappa4 = sel$kappa4[1],
                            parameters = sel$parameters[1])
  }
  names(selections) <- vapply(order_sets, function(x) paste(.mar_check_orders(x), collapse = "+"), character(1))
  summary <- .mar_bind_rows_fill(rows)
  summary <- summary[order(ifelse(is.finite(summary$value), summary$value, Inf)), , drop = FALSE]
  rownames(summary) <- NULL
  list(summary = summary, selections = selections)
}

# ============================================================
# 11. Simulation compatibility helpers
# ============================================================

rskt <- function(n, df = 5, gamma = 1, standardize = TRUE) {
  marpoly_rskewt(n = n, df = df, gamma = gamma, standardize = standardize)
}

simulate_MAR <- function(T, phi_plus = NULL, phi_star = NULL,
                         noise_dist = NULL, noise_args = list(), seed = NULL,
                         burnin = 500, standardize = TRUE,
                         check_stationary = TRUE, ...) {
  T <- .mar_check_positive_integer(T, "T", min_value = 1L)
  phi_plus <- .mar_clean_coefs(phi_plus, "phi_plus")
  phi_star <- .mar_clean_coefs(phi_star, "phi_star")
  if (check_stationary) {
    info <- .mar_stationarity_info(phi_plus, phi_star)
    if (!info$stationary) {
      stop("Stationarity restriction failed: all roots of phi_plus and phi_star must lie outside the unit circle.",
           call. = FALSE)
    }
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit({ if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv) }, add = TRUE)
  }
  innov_fun <- function(n) {
    if (is.null(noise_dist)) {
      stats::rnorm(n)
    } else {
      do.call(noise_dist, c(list(n = n), noise_args))
    }
  }
  sim <- marpoly_simulate(n = T, causal = phi_plus, noncausal = phi_star,
                          innov = innov_fun, burnin = burnin, seed = NULL,
                          standardize = standardize, check_stationary = check_stationary, ...)
  out <- list(X = sim$x, x = sim$x, epsilon = sim$innovations,
              innovations = sim$innovations, phi_plus = phi_plus, phi_star = phi_star,
              T = T, burnin = burnin, stationary = sim$stationary, sim = sim)
  class(out) <- c("simulate_MAR", "list")
  out
}
