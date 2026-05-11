# Estimation and root selection for a MAR(1,2)
library(marpoly)

sim <- marpoly_simulate(
  n = 900,
  causal = 0.25,
  noncausal = c(0.50, 0.20),
  dist = "skewt",
  df = 12,
  gamma = 1.4,
  seed = 6001
)

# The stationary restriction is checked automatically in the simulator.
marpoly_stationary(causal = 0.25, noncausal = c(0.50, 0.20))

fit <- marpoly_fit(
  sim$x,
  p = 3,
  orders = c(2, 3),
  grid_size = 20,
  optimizer = "single",
  maxit = 4000
)

summary(fit)
marpoly_rank(fit, n = 8)
