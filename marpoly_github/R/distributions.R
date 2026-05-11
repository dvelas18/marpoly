#' Draw from the Fernandez-Steel skew Student t distribution
#'
#' The implementation uses the two-piece skewing mechanism. For gamma = 1 the
#' draw reduces to a usual Student t draw. By default, the simulated sample is
#' centered and scaled, which is convenient for MAR simulations.
#'
#' @param n Number of observations.
#' @param df Degrees of freedom.
#' @param gamma Skewness parameter. Values above one produce right skewness;
#'   values below one produce left skewness.
#' @param standardize If TRUE, center and scale the simulated sample.
#' @return Numeric vector.
#' @export
marpoly_rskewt <- function(n, df = 5, gamma = 1, standardize = TRUE) {
  n <- as.integer(n)
  if (n <= 0L) stop("n must be positive.", call. = FALSE)
  if (!is.finite(df) || df <= 0) stop("df must be positive.", call. = FALSE)
  if (!is.finite(gamma) || gamma <= 0) stop("gamma must be positive.", call. = FALSE)
  base <- abs(stats::rt(n, df = df))
  p_pos <- gamma^2 / (1 + gamma^2)
  sign_pos <- stats::runif(n) < p_pos
  x <- numeric(n)
  x[sign_pos] <- gamma * base[sign_pos]
  x[!sign_pos] <- -base[!sign_pos] / gamma
  if (standardize) {
    sx <- stats::sd(x)
    if (is.finite(sx) && sx > 0) {
      x <- (x - mean(x)) / sx
    } else {
      x <- x - mean(x)
    }
  }
  x
}

#' Draw innovation sequences for MAR simulations
#'
#' @param n Number of observations.
#' @param dist Distribution name: "normal", "t", or "skewt".
#' @param df Degrees of freedom for t and skew-t innovations.
#' @param gamma Skewness parameter for skew-t innovations.
#' @param standardize If TRUE, center and scale the generated sample.
#' @param ... Additional arguments reserved for future use.
#' @return Numeric vector of innovations.
#' @export
marpoly_rinnov <- function(n, dist = "normal", df = 5, gamma = 1,
                           standardize = TRUE, ...) {
  if (is.function(dist)) {
    x <- dist(n)
  } else {
    dist <- tolower(as.character(dist)[1L])
    if (dist %in% c("normal", "gaussian", "norm")) {
      x <- stats::rnorm(n)
    } else if (dist %in% c("t", "student", "student-t", "student_t")) {
      x <- stats::rt(n, df = df)
    } else if (dist %in% c("skewt", "skew-t", "skew_t", "fs-skewt")) {
      x <- marpoly_rskewt(n, df = df, gamma = gamma, standardize = FALSE)
    } else {
      stop("Unknown innovation distribution. Use 'normal', 't', 'skewt', or pass a function.", call. = FALSE)
    }
  }
  x <- as.numeric(x)
  if (length(x) != n) {
    stop("The innovation generator must return exactly n observations.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("Innovations must be finite.", call. = FALSE)
  }
  if (standardize) {
    sx <- stats::sd(x)
    if (is.finite(sx) && sx > 0) {
      x <- (x - mean(x)) / sx
    } else {
      x <- x - mean(x)
    }
  }
  x
}
