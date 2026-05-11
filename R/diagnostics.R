#' Residual diagnostics for MAR errors
#'
#' Computes Ljung-Box tests and autocorrelation values for residual powers,
#' commonly the residuals and squared residuals.
#'
#' @param residuals Numeric residual vector.
#' @param lag.max Maximum lag for ACF and Ljung-Box tests.
#' @param powers Integer powers to diagnose.
#' @param fitdf Degrees of freedom subtracted in Ljung-Box tests.
#' @return Object of class `marpoly_diagnostics`.
#' @export
marpoly_diagnostics <- function(residuals, lag.max = 20, powers = c(1, 2), fitdf = 0) {
  e <- as.numeric(residuals)
  e <- e[is.finite(e)]
  if (length(e) < 5L) stop("residuals must contain at least five finite observations.", call. = FALSE)
  lag.max <- min(as.integer(lag.max), length(e) - 1L)
  powers <- as.integer(powers)
  tests <- data.frame(power = powers, statistic = NA_real_, p.value = NA_real_, lag = lag.max)
  acfs <- vector("list", length(powers))
  names(acfs) <- paste0("power", powers)
  for (i in seq_along(powers)) {
    z <- e^powers[i]
    z <- z - mean(z)
    ac <- stats::acf(z, lag.max = lag.max, plot = FALSE, na.action = stats::na.pass)
    acfs[[i]] <- data.frame(lag = as.integer(ac$lag[-1L]), acf = as.numeric(ac$acf[-1L]))
    lb <- stats::Box.test(z, lag = lag.max, type = "Ljung-Box", fitdf = fitdf)
    tests$statistic[i] <- lb$statistic
    tests$p.value[i] <- lb$p.value
  }
  out <- list(tests = tests, acf = acfs, lag.max = lag.max, powers = powers)
  class(out) <- "marpoly_diagnostics"
  out
}

#' @export
print.marpoly_diagnostics <- function(x, ...) {
  cat("marpoly residual diagnostics\n")
  print(x$tests)
  invisible(x)
}
