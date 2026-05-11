# Identification of MAR(1,1) using orders 2 and 4
library(marpoly)

sim <- marpoly_simulate(
  n = 800,
  causal = 0.20,
  noncausal = 0.70,
  dist = "t",
  df = 5,
  seed = 1002
)

fit <- marpoly_fit(
  sim$x,
  p = 2,
  orders = c(2, 4),
  grid_size = 12,
  optimizer = "single"
)

summary(fit)
marpoly_rank(fit)
plot(fit, which = "rank")
