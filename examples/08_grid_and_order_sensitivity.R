# Sensitivity to frequency grid size and spectral order combinations
library(marpoly)

sim <- marpoly_simulate(
  n = 700,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 8,
  gamma = 1.5,
  seed = 8001
)

order_sets <- list(c(2, 3), c(2, 4), c(2, 3, 4))
fits <- lapply(order_sets, function(ord) {
  marpoly_fit(sim$x, p = 2, orders = ord, grid_size = if (4 %in% ord) 12 else 24)
})

lapply(fits, marpoly_rank, n = 4)
