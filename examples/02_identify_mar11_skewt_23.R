# Identification of MAR(1,1) using orders 2 and 3
library(marpoly)

sim <- marpoly_simulate(
  n = 800,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 20,
  gamma = 1.5,
  seed = 1001
)

fit <- marpoly_fit(
  sim$x,
  p = 2,
  orders = c(2, 3),
  grid_size = 24,
  optimizer = "single"
)

summary(fit)
marpoly_rank(fit)
plot(fit, which = "diagnostics")
