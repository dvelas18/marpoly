# Errors for each root allocation and ACF checks
library(marpoly)

sim <- simulate_MAR(
  T = 5000,
  phi_plus = 0.20,
  phi_star = 0.70,
  noise_dist = rskt,
  noise_args = list(df = 20, gamma = 1.5),
  seed = 7001
)

fit <- marpoly_fit(
  sim$X,
  p = 2,
  orders = c(2, 3),
  control = list(
    max_bispec_pairs = 50000,
    n_random_starts = 5,
    compute_errors = TRUE,
    compute_se = TRUE
  )
)

summary(fit)
marpoly_rank(fit)

# Selected allocation.
e_selected <- marpoly_errors(fit)
marpoly_diagnostics(e_selected, powers = c(1, 2), lag.max = 20)

# Every root allocation.
e_all <- marpoly_errors(fit, allocation = "all")
diag_all <- lapply(e_all, marpoly_diagnostics, powers = c(1, 2), lag.max = 20)
diag_all

plot(fit, which = "diagnostics")
