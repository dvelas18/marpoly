#' Quick plots for series, periodograms, spectra, and fitted MAR models
#'
#' @param x Numeric time series.
#' @param orders Spectral orders to display.
#' @param fit Optional `marpoly_fit` object. When supplied, the order-2 fitted
#'   spectral density is overlaid on the periodogram.
#' @param grid_size Frequency grid size.
#' @param component Component used for order 3 and 4 heatmaps.
#' @param slice Slice used for order 4 heatmaps.
#' @param demean If TRUE, demean the series before computing periodograms.
#' @return Invisibly returns the computed periodograms.
#' @export
marpoly_plot <- function(x, orders = c(2, 3), fit = NULL, grid_size = NULL,
                         component = c("real", "imaginary", "modulus"),
                         slice = 1L, demean = TRUE) {
  component <- match.arg(component)
  x_clean <- mp_clean_series(x, demean = demean)
  orders <- mp_parse_orders(orders)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  n_panels <- 1L + length(orders)
  graphics::par(mfrow = c(n_panels, 1L), mar = c(4, 4, 3, 1))
  graphics::plot(seq_along(x_clean), x_clean, type = "l", xlab = "Time", ylab = "Series", main = "Time series")
  pgs <- list()
  for (k in orders) {
    pg <- marpoly_periodogram(x_clean, order = k, grid_size = grid_size, demean = FALSE)
    pgs[[paste0("order", k)]] <- pg
    if (k == 2L) {
      graphics::plot(pg$frequency, Re(pg$values), type = "l", xlab = "Frequency", ylab = "I2", main = "Periodogram, order 2")
      if (!is.null(fit) && (inherits(fit, "marpoly_fit") || inherits(fit, "MAR_selection"))) {
        sf <- mp_selected_fit(fit)
        kap2 <- if (!is.null(sf$kappa2) && is.finite(sf$kappa2)) sf$kappa2 else 1
        sp <- marpoly_spectrum(causal = sf$phi_plus, noncausal = sf$phi_star, order = 2, kappa = kap2, omega = pg$omega)
        graphics::lines(pg$frequency, Re(sp$values), lwd = 2)
        graphics::legend("topright", legend = c("periodogram", "fitted spectrum"), lty = c(1, 1), bty = "n")
      }
    } else {
      plot(pg, component = component, slice = slice)
    }
  }
  invisible(pgs)
}

#' @export
plot.marpoly_fit <- function(x, which = c("rank", "diagnostics", "periodogram"), ...) {
  which <- match.arg(which)
  if (which == "rank") {
    nn <- if (!is.null(x$all_root_allocations)) nrow(x$all_root_allocations) else if (!is.null(x$results)) nrow(x$results) else 10L
    tab <- marpoly_rank(x, n = nn, root_allocations = TRUE)
    graphics::plot(tab$rank, tab$objective, type = "b", xlab = "Rank", ylab = "Objective", main = "Ranked root allocations")
  } else if (which == "diagnostics") {
    e <- marpoly_errors(x)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = c(2, 1))
    stats::acf(e, main = "ACF of residuals", ...)
    stats::acf(e^2, main = "ACF of squared residuals", ...)
  } else {
    marpoly_plot(x$x, orders = x$orders %||% mp_selected_fit(x)$orders, fit = x, demean = FALSE, ...)
  }
  invisible(x)
}
