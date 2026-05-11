# Estimation and root selection for a MAR(1,1)
library(marpoly)

sim <- simulate_MAR(
  T = 5000,
  phi_plus = 0.20,
  phi_star = 0.70,
  noise_dist = rskt,
  noise_args = list(df = 4, gamma = 1.5),
  seed = 5001
)

fit_11 <- MAR_spectral_estimate(
  y = sim$X,
  r = 1,
  s = 1,
  orders = c(2, 3, 4),
  control = list(
    max_bispec_pairs = 50000,
    max_trispec_tuples = 10000,
    trispec_max_index = 35,
    n_random_starts = 5,
    compute_errors = TRUE,
    compute_se = TRUE
  )
)

fit_11$parameter_table
fit_11$value
fit_11$R2
fit_11$R3
fit_11$R4
fit_11$error_diagnostics

selection <- MAR_select_all_by_second_order(
  y = sim$X,
  orders = c(2, 3, 4),
  max_p = 5,
  ic = "bic",
  control = list(
    max_bispec_pairs = 50000,
    max_trispec_tuples = 10000,
    trispec_max_index = 35,
    n_random_starts = 5,
    compute_errors = TRUE,
    compute_se = TRUE
  )
)

selection$selected_p
selection$eta_pseudo_causal
selection$report
selection$selected
selection$all_root_allocations
