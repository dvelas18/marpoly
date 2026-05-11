# Estimation and root selection for a MAR(1,2)
library(marpoly)

sim <- simulate_MAR(
  T = 5000,
  phi_plus = 0.25,
  phi_star = c(0.50, 0.20),
  noise_dist = rskt,
  noise_args = list(df = 12, gamma = 1.4),
  seed = 6001
)

selection <- MAR_select_all_by_second_order(
  y = sim$X,
  orders = c(2, 3),
  max_p = 5,
  ic = "bic",
  control = list(
    max_bispec_pairs = 50000,
    n_random_starts = 5,
    compute_errors = TRUE,
    compute_se = TRUE
  )
)

selection$report
selection$selected
selection$all_root_allocations
