#' OLS pseudo-causal AR(p) estimation
#'
#' Fits the pseudo-causal AR(p) representation by ordinary least squares. This
#' is used for second-order order selection and for initializing root allocations.
#'
#' @param x Numeric time series.
#' @param p AR order.
#' @param demean If TRUE, subtract the sample mean before fitting.
#' @param include_mean If TRUE, include an intercept in the AR regression.
#' @return Object of class `marpoly_ar`.
#' @export
marpoly_ar_ols <- function(x, p, demean = TRUE, include_mean = FALSE) {
  x <- mp_clean_series(x, demean = demean)
  p <- as.integer(p)
  if (p < 1L) stop("p must be at least one.", call. = FALSE)
  n <- length(x)
  if (n <= p + 2L) stop("Not enough observations for the requested AR order.", call. = FALSE)
  emb <- stats::embed(x, p + 1L)
  y <- emb[, 1L]
  X <- emb[, -1L, drop = FALSE]
  colnames(X) <- paste0("ar", seq_len(p))
  dat <- data.frame(y = y, X)
  if (include_mean) {
    fit <- stats::lm(y ~ ., data = dat)
  } else {
    fit <- stats::lm(y ~ . - 1, data = dat)
  }
  cf <- stats::coef(fit)
  ar_names <- paste0("ar", seq_len(p))
  ar <- as.numeric(cf[ar_names])
  names(ar) <- ar_names
  res <- stats::residuals(fit)
  sigma2 <- mean(res^2)
  k <- length(cf)
  loglik <- -0.5 * length(res) * (log(2 * pi) + log(sigma2) + 1)
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(length(res)) * k
  out <- list(
    ar = ar,
    coefficients = cf,
    fit = fit,
    residuals = res,
    fitted = stats::fitted(fit),
    p = p,
    nobs = length(res),
    sigma2 = sigma2,
    sigma = sqrt(sigma2),
    logLik = loglik,
    AIC = aic,
    BIC = bic,
    demean = demean,
    include_mean = include_mean,
    roots = mp_roots(ar)
  )
  class(out) <- "marpoly_ar"
  out
}

#' Select pseudo-causal AR order by information criteria
#'
#' @param x Numeric time series.
#' @param p_max Maximum AR order.
#' @param criterion Selection criterion: "BIC" or "AIC".
#' @param demean If TRUE, subtract the sample mean before fitting.
#' @param include_mean If TRUE, include an intercept in each AR regression.
#' @return Object of class `marpoly_ar_selection`.
#' @export
marpoly_select_p <- function(x, p_max = 6, criterion = c("BIC", "AIC"),
                             demean = TRUE, include_mean = FALSE) {
  criterion <- match.arg(criterion)
  p_max <- as.integer(p_max)
  if (p_max < 1L) stop("p_max must be at least one.", call. = FALSE)
  fits <- vector("list", p_max)
  table <- data.frame(p = seq_len(p_max), AIC = NA_real_, BIC = NA_real_, sigma = NA_real_,
                      LjungBox_pvalue = NA_real_)
  for (p in seq_len(p_max)) {
    fits[[p]] <- marpoly_ar_ols(x, p = p, demean = demean, include_mean = include_mean)
    table$AIC[p] <- fits[[p]]$AIC
    table$BIC[p] <- fits[[p]]$BIC
    table$sigma[p] <- fits[[p]]$sigma
    lb <- try(stats::Box.test(fits[[p]]$residuals, lag = min(10L, max(1L, length(fits[[p]]$residuals) - 1L)),
                              type = "Ljung-Box", fitdf = p), silent = TRUE)
    if (!inherits(lb, "try-error")) table$LjungBox_pvalue[p] <- lb$p.value
  }
  selected_p <- table$p[which.min(table[[criterion]])]
  out <- list(
    selected_p = selected_p,
    criterion = criterion,
    table = table,
    fits = fits,
    selected_fit = fits[[selected_p]]
  )
  class(out) <- "marpoly_ar_selection"
  out
}

#' @export
coef.marpoly_ar <- function(object, ...) {
  object$coefficients
}

#' @export
residuals.marpoly_ar <- function(object, ...) {
  object$residuals
}

#' @export
fitted.marpoly_ar <- function(object, ...) {
  object$fitted
}

#' @export
summary.marpoly_ar <- function(object, ...) {
  sm <- summary(object$fit)
  tab <- sm$coefficients
  out <- list(
    coefficients = tab,
    sigma = object$sigma,
    nobs = object$nobs,
    p = object$p,
    AIC = object$AIC,
    BIC = object$BIC,
    roots = object$roots
  )
  class(out) <- "summary_marpoly_ar"
  out
}

#' @export
print.marpoly_ar <- function(x, ...) {
  cat("Pseudo-causal AR(", x$p, ") by OLS\n", sep = "")
  print(summary(x)$coefficients)
  cat("sigma:", round(x$sigma, 6), " AIC:", round(x$AIC, 3), " BIC:", round(x$BIC, 3), "\n")
  invisible(x)
}

#' @export
print.summary_marpoly_ar <- function(x, ...) {
  cat("Pseudo-causal AR(", x$p, ") regression summary\n", sep = "")
  print(x$coefficients)
  cat("sigma:", round(x$sigma, 6), " nobs:", x$nobs, " AIC:", round(x$AIC, 3), " BIC:", round(x$BIC, 3), "\n")
  invisible(x)
}

#' @export
print.marpoly_ar_selection <- function(x, ...) {
  cat("Pseudo-causal AR order selection\n")
  cat("  criterion:  ", x$criterion, "\n", sep = "")
  cat("  selected p: ", x$selected_p, "\n", sep = "")
  print(x$table)
  invisible(x)
}
