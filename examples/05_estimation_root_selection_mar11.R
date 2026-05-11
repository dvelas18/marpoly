# Estimation and root selection for a MAR(1,1)
library(marpoly)

sim <- marpoly_simulate(
  n = 700,
  causal = 0.25,
  noncausal = 0.75,
  dist = "skewt",
  df = 10,
  gamma = 1.5,
  seed = 5001
)

# Second-order order selection.
sel <- marpoly_select_p(sim$x, p_max = 5, criterion = "BIC")
print(sel)
summary(sel$selected_fit)

# Estimate all possible allocations for the selected p.
fit <- marpoly_fit(sim$x, p = sel$selected_p, orders = c(2, 3), grid_size = 24)
summary(fit)

# The compact ranking shows the selected model and main competitors.
marpoly_rank(fit)
