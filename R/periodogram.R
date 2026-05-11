#' Discrete Fourier transform for a time series
#'
#' @param x Numeric time series.
#' @param demean If TRUE, subtract the sample mean before applying the DFT.
#' @return A list with DFT values, angular frequencies, and unit frequencies.
#' @export
marpoly_dft <- function(x, demean = TRUE) {
  x <- mp_clean_series(x, demean = demean)
  n <- length(x)
  d <- .mar_dft(x)
  list(
    dft = d,
    omega = 2 * pi * (0:(n - 1L)) / n,
    frequency = (0:(n - 1L)) / n,
    n = n,
    demean = demean
  )
}

#' Higher-order periodogram
#'
#' Computes the periodogram of order 2, the biperiodogram of order 3, or a
#' fourth-order periodogram used for trispectral estimation.
#'
#' @param x Numeric time series.
#' @param order Periodogram order: 2, 3, or 4.
#' @param grid_size Number of frequencies used in each frequency dimension.
#' @param demean If TRUE, subtract the sample mean.
#' @param indices Optional integer Fourier indices. If supplied, `grid_size` is
#'   ignored.
#' @return Object of class `marpoly_periodogram`.
#' @export
marpoly_periodogram <- function(x, order = 2, grid_size = NULL, demean = TRUE,
                                indices = NULL) {
  order <- mp_order_check(order)
  x <- mp_clean_series(x, demean = demean)
  n <- length(x)
  d <- .mar_dft(x)
  if (is.null(indices)) {
    indices <- mp_frequency_indices(n, grid_size = grid_size, order = order, positive = order == 2L)
  } else {
    indices <- as.integer(indices)
    indices <- indices[indices >= 1L & indices <= n - 1L]
    if (length(indices) == 0L) stop("indices must include at least one Fourier index between 1 and n - 1.", call. = FALSE)
  }
  m <- length(indices)
  omega <- 2 * pi * indices / n
  freq <- indices / n
  if (order == 2L) {
    values <- d[indices + 1L] * d[mp_mod(-indices, n) + 1L] / (2 * pi * n)
    values <- Re(values)
  } else if (order == 3L) {
    values <- matrix(0+0i, m, m)
    for (a in seq_len(m)) {
      for (b in seq_len(m)) {
        k <- mp_mod(-indices[a] - indices[b], n)
        values[a, b] <- d[indices[a] + 1L] * d[indices[b] + 1L] * d[k + 1L] / ((2 * pi)^2 * n)
      }
    }
    dimnames(values) <- list(freq, freq)
  } else {
    values <- array(0+0i, dim = c(m, m, m), dimnames = list(freq, freq, freq))
    for (a in seq_len(m)) {
      for (b in seq_len(m)) {
        for (c in seq_len(m)) {
          k <- mp_mod(-indices[a] - indices[b] - indices[c], n)
          values[a, b, c] <- d[indices[a] + 1L] * d[indices[b] + 1L] * d[indices[c] + 1L] * d[k + 1L] / ((2 * pi)^3 * n)
        }
      }
    }
  }
  out <- list(
    values = values,
    frequency = freq,
    omega = omega,
    indices = indices,
    order = order,
    n = n,
    demean = demean,
    series = x
  )
  class(out) <- "marpoly_periodogram"
  out
}

#' Compute multiple periodograms
#'
#' @param x Numeric time series.
#' @param orders Vector containing any subset of 2, 3, and 4.
#' @param ... Arguments passed to `marpoly_periodogram`.
#' @return Named list of periodogram objects.
#' @export
marpoly_periodograms <- function(x, orders = c(2, 3, 4), ...) {
  orders <- mp_parse_orders(orders)
  out <- lapply(orders, function(k) marpoly_periodogram(x, order = k, ...))
  names(out) <- paste0("order", orders)
  out
}

#' @export
plot.marpoly_periodogram <- function(x, component = c("real", "imaginary", "modulus"),
                                     slice = 1L, main = NULL, ...) {
  component <- match.arg(component)
  val <- x$values
  get_component <- function(z) {
    if (component == "real") Re(z) else if (component == "imaginary") Im(z) else Mod(z)
  }
  if (is.null(main)) {
    main <- paste0("Order ", x$order, " periodogram")
  }
  if (x$order == 2L) {
    graphics::plot(x$frequency, Re(val), type = "l", xlab = "Frequency", ylab = "I2", main = main, ...)
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
