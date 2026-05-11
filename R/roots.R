#' Check MAR stationarity restrictions
#'
#' For the MAR(r,s) convention used in this package, both the causal lag
#' polynomial and the noncausal lead polynomial must have roots outside the
#' unit circle.
#'
#' @param causal Numeric vector of causal lag coefficients.
#' @param noncausal Numeric vector of noncausal lead coefficients.
#' @param tol Numerical tolerance used in the unit-root boundary check.
#' @return A list with logical result and polynomial roots.
#' @export
marpoly_stationary <- function(causal = NULL, noncausal = NULL, tol = 1e-8) {
  causal <- mp_as_vector(causal, "causal")
  noncausal <- mp_as_vector(noncausal, "noncausal")
  roots_causal <- if (length(causal)) mp_roots(causal) else complex(0)
  roots_noncausal <- if (length(noncausal)) mp_roots(noncausal) else complex(0)
  ok_causal <- length(causal) == 0L || all(Mod(roots_causal) > 1 + tol)
  ok_noncausal <- length(noncausal) == 0L || all(Mod(roots_noncausal) > 1 + tol)
  out <- list(
    stationary = ok_causal && ok_noncausal,
    causal_ok = ok_causal,
    noncausal_ok = ok_noncausal,
    causal_roots = roots_causal,
    noncausal_roots = roots_noncausal,
    tol = tol
  )
  class(out) <- "marpoly_stationary"
  out
}

mp_check_stationary <- function(causal = NULL, noncausal = NULL, tol = 1e-8) {
  st <- marpoly_stationary(causal, noncausal, tol = tol)
  if (!st$stationary) {
    msg <- "Stationarity restriction failed: roots of both MAR polynomials must lie outside the unit circle."
    if (!st$causal_ok) {
      msg <- paste0(msg, " Causal roots: ", paste(round(Mod(st$causal_roots), 4), collapse = ", "), ".")
    }
    if (!st$noncausal_ok) {
      msg <- paste0(msg, " Noncausal roots: ", paste(round(Mod(st$noncausal_roots), 4), collapse = ", "), ".")
    }
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

mp_root_groups <- function(lambda, tol = 1e-7) {
  lambda <- as.complex(lambda)
  p <- length(lambda)
  used <- rep(FALSE, p)
  groups <- list()
  if (p == 0L) {
    return(groups)
  }
  for (i in seq_len(p)) {
    if (used[i]) next
    if (abs(Im(lambda[i])) < tol) {
      groups[[length(groups) + 1L]] <- i
      used[i] <- TRUE
    } else {
      candidates <- which(!used & seq_len(p) != i & abs(lambda - Conj(lambda[i])) < tol)
      if (length(candidates) > 0L) {
        j <- candidates[1L]
        groups[[length(groups) + 1L]] <- c(i, j)
        used[c(i, j)] <- TRUE
      } else {
        groups[[length(groups) + 1L]] <- i
        used[i] <- TRUE
      }
    }
  }
  groups
}

#' Enumerate root allocations for a pseudo-causal AR(p)
#'
#' This function factorizes the pseudo-causal AR polynomial and enumerates all
#' feasible allocations of inverse roots between the causal and noncausal
#' polynomials. Complex conjugate roots are allocated together so that real
#' model coefficients are preserved whenever possible.
#'
#' @param ar_coef Numeric vector of pseudo-causal AR coefficients.
#' @param tol Numerical tolerance for real/conjugate detection.
#' @return A data frame with one row per allocation and list-columns containing
#'   coefficients and inverse roots.
#' @export
marpoly_root_allocations <- function(ar_coef, tol = 1e-7) {
  ar_coef <- mp_as_vector(ar_coef, "ar_coef")
  p <- length(ar_coef)
  if (p == 0L) {
    stop("ar_coef must contain at least one AR coefficient.", call. = FALSE)
  }
  roots <- mp_roots(ar_coef)
  lambda <- 1 / roots
  groups <- mp_root_groups(lambda, tol = tol)
  masks <- mp_all_binary(length(groups))
  rows <- vector("list", nrow(masks))
  for (i in seq_len(nrow(masks))) {
    causal_idx <- integer(0)
    noncausal_idx <- integer(0)
    for (g in seq_along(groups)) {
      if (masks[i, g] == 1L) {
        causal_idx <- c(causal_idx, groups[[g]])
      } else {
        noncausal_idx <- c(noncausal_idx, groups[[g]])
      }
    }
    causal_lam <- lambda[causal_idx]
    noncausal_lam <- lambda[noncausal_idx]
    causal <- mp_coefs_from_lambdas(causal_lam, tol = tol)
    noncausal <- mp_coefs_from_lambdas(noncausal_lam, tol = tol)
    r <- length(causal)
    s <- length(noncausal)
    rows[[i]] <- data.frame(
      allocation_id = i,
      model = mp_model_label(r, s),
      r = r,
      s = s,
      stringsAsFactors = FALSE
    )
    rows[[i]]$causal <- list(as.numeric(Re(causal)))
    rows[[i]]$noncausal <- list(as.numeric(Re(noncausal)))
    rows[[i]]$causal_inverse_roots <- list(causal_lam)
    rows[[i]]$noncausal_inverse_roots <- list(noncausal_lam)
  }
  out <- do.call(rbind, rows)
  out$order_p <- p
  out
}
