# Periodograms, spectral densities, and plots
library(marpoly)

sim <- marpoly_simulate(
  n = 500,
  causal = 0.25,
  noncausal = 0.65,
  dist = "skewt",
  df = 8,
  gamma = 1.4,
  seed = 42
)

# Empirical periodograms of orders 2, 3, and 4.
p2 <- marpoly_periodogram(sim$x, order = 2)
p3 <- marpoly_periodogram(sim$x, order = 3, grid_size = 32)
p4 <- marpoly_periodogram(sim$x, order = 4, grid_size = 12)

plot(p2)
plot(p3, component = "real")
plot(p3, component = "imaginary")
plot(p4, component = "modulus", slice = 4)

# Theoretical spectra under the DGP parameters.
s2 <- marpoly_spectrum(causal = 0.25, noncausal = 0.65, order = 2, kappa = 1, omega = p2$omega)
s3 <- marpoly_spectrum(causal = 0.25, noncausal = 0.65, order = 3, kappa = 1, grid_size = 32)

plot(s2)
plot(s3, component = "real")

# Combined plotting helper.
marpoly_plot(sim$x, orders = c(2, 3), grid_size = 32)
