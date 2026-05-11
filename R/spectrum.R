#' MAR transfer function
#'
#' Computes the transfer function for
#' `(1 - causal(L))(1 - noncausal(L^-1)) y_t = epsilon_t`.
#'
#' @param omega Angular frequencies in radians.
#' @param causal Numeric vector of causal lag coefficients.
#' @param noncausal Numeric vector of noncausal lead coefficients.
#' @return Complex-valued transfer function.
#' @export
marpoly_transfer <- function(omega, causal = NULL, noncausal = NULL) {
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  z <- exp(-1i * omega)
  den_c <- 1 + 0i + 0 * z
  den_n <- 1 + 0i + 0 * z
  if (length(causal) > 0L) {
    for (j in seq_along(causal)) {
      den_c <- den_c - causal[j] * z^j
    }
  }
  if (length(noncausal) > 0L) {
    for (j in seq_along(noncausal)) {
      den_n <- den_n - noncausal[j] * z^(-j)
    }
  }
  1 / (den_c * den_n)
}

mp_transfer_product <- function(omega, order, causal = NULL, noncausal = NULL) {
  order <- mp_order_check(order)
  w <- omega
  psi_w <- marpoly_transfer(w, causal = causal, noncausal = noncausal)
  if (order == 2L) {
    prod <- psi_w * marpoly_transfer(-w, causal = causal, noncausal = noncausal)
  } else if (order == 3L) {
    base <- outer(psi_w, psi_w, "*")
    prod <- base * marpoly_transfer(-outer(w, w, "+"), causal = causal, noncausal = noncausal)
  } else {
    m <- length(w)
    prod <- array(0+0i, dim = c(m, m, m))
    base <- outer(psi_w, psi_w, "*")
    sum12 <- outer(w, w, "+")
    for (k in seq_len(m)) {
      prod[, , k] <- base * psi_w[k] * marpoly_transfer(-(sum12 + w[k]), causal = causal, noncausal = noncausal)
    }
  }
  prod
}

#' Theoretical MAR polyspectral density
#'
#' Computes the spectrum (order 2), bispectrum (order 3), or trispectrum-type
#' fourth-order spectral density implied by a MAR(r,s) model and an innovation
#' cumulant.
#'
#' @param causal Numeric vector of causal lag coefficients.
#' @param noncausal Numeric vector of noncausal lead coefficients.
#' @param order Spectral order: 2, 3, or 4.
#' @param kappa Innovation cumulant of the selected order.
#' @param grid_size Number of frequencies per dimension when `omega` is NULL.
#' @param omega Optional angular frequency grid.
#' @param n Frequency denominator used only to construct a grid when `omega` is
#'   NULL.
#' @return Object of class `marpoly_spectrum`.
#' @export
marpoly_spectrum <- function(causal = NULL, noncausal = NULL, order = 2,
                             kappa = 1, grid_size = NULL, omega = NULL,
                             n = 512) {
  order <- mp_order_check(order)
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  if (is.null(omega)) {
    idx <- mp_frequency_indices(n, grid_size = grid_size, order = order, positive = order == 2L)
    omega <- 2 * pi * idx / n
    freq <- idx / n
  } else {
    omega <- as.numeric(omega)
    freq <- omega / (2 * pi)
  }
  prod <- mp_transfer_product(omega, order = order, causal = causal, noncausal = noncausal)
  values <- kappa / (2 * pi)^(order - 1L) * prod
  if (order == 2L) values <- Re(values)
  out <- list(
    values = values,
    frequency = freq,
    omega = omega,
    order = order,
    kappa = kappa,
    causal = causal,
    noncausal = noncausal,
    model = mp_model_label(length(causal), length(noncausal))
  )
  class(out) <- "marpoly_spectrum"
  out
}

#' Compute multiple theoretical polyspectra
#'
#' @param causal Numeric vector of causal coefficients.
#' @param noncausal Numeric vector of noncausal coefficients.
#' @param orders Vector containing any subset of 2, 3, and 4.
#' @param kappa Optional named vector or list of cumulants.
#' @param ... Additional arguments passed to `marpoly_spectrum`.
#' @return Named list of spectral-density objects.
#' @export
marpoly_spectra <- function(causal = NULL, noncausal = NULL, orders = c(2, 3, 4),
                            kappa = NULL, ...) {
  orders <- mp_parse_orders(orders)
  if (is.null(kappa)) {
    kappa <- rep(1, length(orders))
    names(kappa) <- as.character(orders)
  }
  out <- lapply(orders, function(k) {
    kk <- if (!is.null(names(kappa)) && as.character(k) %in% names(kappa)) kappa[[as.character(k)]] else kappa[[which(orders == k)[1L]]]
    marpoly_spectrum(causal = causal, noncausal = noncausal, order = k, kappa = kk, ...)
  })
  names(out) <- paste0("order", orders)
  out
}

#' @export
plot.marpoly_spectrum <- function(x, component = c("real", "imaginary", "modulus"),
                                  slice = 1L, main = NULL, ...) {
  component <- match.arg(component)
  val <- x$values
  get_component <- function(z) {
    if (component == "real") Re(z) else if (component == "imaginary") Im(z) else Mod(z)
  }
  if (is.null(main)) {
    main <- paste0(x$model, " order ", x$order, " spectral density")
  }
  if (x$order == 2L) {
    graphics::plot(x$frequency, Re(val), type = "l", xlab = "Frequency", ylab = "S2", main = main, ...)
  } else if (x$order == 3L) {
    graphics::image(x$frequency, x$frequency, get_component(val), xlab = expression(omega[1]),
                    ylab = expression(omega[2]), main = paste(main, component), ...)
  } else {
    slice <- max(1L, min(as.integer(slice), dim(val)[3L]))
    graphics::image(x$frequency, x$frequency, get_component(val[, , slice]),
                    xlab = expression(omega[1]), ylab = expression(omega[2]),
                    main = paste(main, component, "slice", slice), ...)
  }
  invisible(x)
}
