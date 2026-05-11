# Errors for each root allocation and ACF checks
library(marpoly)

sim <- marpoly_simulate(
  n = 700,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 20,
  gamma = 1.5,
  seed = 7001
)

fit <- marpoly_fit(sim$x, p = 2, orders = c(2, 3), grid_size = 24)
summary(fit)

# Selected allocation.
e_selected <- marpoly_errors(fit)
marpoly_diagnostics(e_selected, powers = c(1, 2), lag.max = 20)

# Every allocation.
e_all <- marpoly_errors(fit, allocation = "all")
diag_all <- lapply(e_all, marpoly_diagnostics, powers = c(1, 2), lag.max = 20)
diag_all

plot(fit, which = "diagnostics")
