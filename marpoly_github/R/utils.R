# Internal utility functions for marpoly

mp_as_vector <- function(x, name = "x") {
  if (is.null(x)) {
    return(numeric(0))
  }
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite numeric values.", call. = FALSE)
  }
  x
}

mp_clean_series <- function(x, demean = TRUE) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 5L) {
    stop("x must contain at least five finite observations.", call. = FALSE)
  }
  if (demean) {
    x <- x - mean(x)
  }
  x
}

mp_mod <- function(a, n) {
  a %% n
}

mp_poly_multiply <- function(a, b) {
  a <- as.complex(a)
  b <- as.complex(b)
  out <- rep(0+0i, length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    for (j in seq_along(b)) {
      out[i + j - 1L] <- out[i + j - 1L] + a[i] * b[j]
    }
  }
  out
}

mp_real_if_close <- function(z, tol = 1e-8) {
  if (length(z) == 0L) {
    return(numeric(0))
  }
  if (all(abs(Im(z)) < tol)) {
    return(as.numeric(Re(z)))
  }
  z
}

mp_coefs_from_lambdas <- function(lambda, tol = 1e-8) {
  lambda <- as.complex(lambda)
  if (length(lambda) == 0L) {
    return(numeric(0))
  }
  poly <- 1+0i
  for (lam in lambda) {
    poly <- mp_poly_multiply(poly, c(1+0i, -lam))
  }
  coefs <- -poly[-1L]
  mp_real_if_close(coefs, tol = tol)
}

mp_roots <- function(coefs) {
  coefs <- mp_as_vector(coefs, "coefs")
  if (length(coefs) == 0L) {
    return(complex(0))
  }
  stats::polyroot(c(1, -coefs))
}

mp_is_stationary_poly <- function(coefs, tol = 1e-8) {
  coefs <- mp_as_vector(coefs, "coefs")
  if (length(coefs) == 0L) {
    return(TRUE)
  }
  roots <- mp_roots(coefs)
  all(is.finite(Mod(roots))) && all(Mod(roots) > 1 + tol)
}

mp_model_label <- function(r, s) {
  if (r > 0L && s == 0L) {
    paste0("AR(", r, ",0)")
  } else if (r == 0L && s > 0L) {
    paste0("AR(0,", s, ")")
  } else if (r > 0L && s > 0L) {
    paste0("MAR(", r, ",", s, ")")
  } else {
    "white-noise"
  }
}

mp_format_vector <- function(x, digits = 4L) {
  if (length(x) == 0L) {
    return("")
  }
  paste(format(round(as.numeric(Re(x)), digits), nsmall = min(2L, digits)), collapse = ", ")
}

mp_all_binary <- function(n) {
  if (n == 0L) {
    return(matrix(integer(0), nrow = 1L))
  }
  out <- matrix(0L, nrow = 2^n, ncol = n)
  for (i in seq_len(2^n)) {
    bits <- as.integer(intToBits(i - 1L))[seq_len(n)]
    out[i, ] <- bits
  }
  out
}

mp_numeric_hessian <- function(fn, x, eps = 1e-4) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) {
    return(matrix(numeric(0), 0L, 0L))
  }
  H <- matrix(0, n, n)
  fx <- fn(x)
  step <- eps * pmax(abs(x), 1)
  for (i in seq_len(n)) {
    ei <- rep(0, n)
    ei[i] <- step[i]
    fpi <- fn(x + ei)
    fmi <- fn(x - ei)
    H[i, i] <- (fpi - 2 * fx + fmi) / (step[i]^2)
    if (i < n) {
      for (j in (i + 1L):n) {
        ej <- rep(0, n)
        ej[j] <- step[j]
        fpp <- fn(x + ei + ej)
        fpm <- fn(x + ei - ej)
        fmp <- fn(x - ei + ej)
        fmm <- fn(x - ei - ej)
        H[i, j] <- (fpp - fpm - fmp + fmm) / (4 * step[i] * step[j])
        H[j, i] <- H[i, j]
      }
    }
  }
  H
}

mp_try_se <- function(fn, par, nobs) {
  out <- rep(NA_real_, length(par))
  if (length(par) == 0L) {
    return(out)
  }
  H <- try(mp_numeric_hessian(fn, par), silent = TRUE)
  if (inherits(H, "try-error") || any(!is.finite(H))) {
    return(out)
  }
  inv <- try(solve(H), silent = TRUE)
  if (inherits(inv, "try-error") || any(!is.finite(inv))) {
    inv <- try(solve(H + diag(1e-6, nrow(H))), silent = TRUE)
  }
  if (inherits(inv, "try-error") || any(!is.finite(inv))) {
    return(out)
  }
  val <- diag(inv) / max(nobs, 1L)
  val[val < 0] <- NA_real_
  sqrt(val)
}

mp_default_grid <- function(n, order) {
  if (order == 2L) {
    min(128L, max(8L, floor((n - 1L) / 2L)))
  } else if (order == 3L) {
    min(48L, max(8L, n - 1L))
  } else {
    min(18L, max(6L, n - 1L))
  }
}

mp_frequency_indices <- function(n, grid_size = NULL, order = 2L, positive = FALSE) {
  if (is.null(grid_size)) {
    grid_size <- mp_default_grid(n, order)
  }
  if (positive) {
    pool <- seq_len(max(1L, floor((n - 1L) / 2L)))
  } else {
    pool <- seq_len(n - 1L)
  }
  grid_size <- min(as.integer(grid_size), length(pool))
  unique(as.integer(round(seq(min(pool), max(pool), length.out = grid_size))))
}

mp_order_check <- function(order) {
  order <- as.integer(order)
  if (length(order) != 1L || !(order %in% c(2L, 3L, 4L))) {
    stop("order must be one of 2, 3, or 4.", call. = FALSE)
  }
  order
}

mp_parse_orders <- function(orders) {
  orders <- sort(unique(as.integer(orders)))
  if (!all(orders %in% c(2L, 3L, 4L))) {
    stop("orders must contain only 2, 3, and/or 4.", call. = FALSE)
  }
  orders
}

mp_ar_filter_forward <- function(e, coefs) {
  coefs <- mp_as_vector(coefs, "coefs")
  e <- as.numeric(e)
  n <- length(e)
  p <- length(coefs)
  y <- numeric(n)
  if (p == 0L) {
    return(e)
  }
  for (t in seq_len(n)) {
    val <- e[t]
    for (j in seq_len(p)) {
      if (t - j >= 1L) {
        val <- val + coefs[j] * y[t - j]
      }
    }
    y[t] <- val
  }
  y
}
