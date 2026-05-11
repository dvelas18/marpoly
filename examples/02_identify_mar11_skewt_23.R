# Identification of MAR(1,1) using orders 2 and 3
library(marpoly)

sim <- simulate_MAR(
  T = 5000,
  phi_plus = 0.20,
  phi_star = 0.70,
  noise_dist = rskt,
  noise_args = list(df = 20, gamma = 1.5),
  seed = 1001
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

selection$selected_p
selection$report
selection$selected
selected_name <- selection$selected$model[1]
selection$fits[[selected_name]]$parameter_table
