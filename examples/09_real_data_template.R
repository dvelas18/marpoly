# Real-data template
# Replace x with a stationary transformation of the empirical series.
library(marpoly)

# Example placeholder: x <- diff(log(price_series))
# x <- x[is.finite(x)]

# sel <- marpoly_select_p(x, p_max = 6, criterion = "BIC")
# print(sel)
#
# fit <- marpoly_fit(x, p = sel$selected_p, orders = c(2, 3, 4), grid_size = 12)
# summary(fit)
# marpoly_rank(fit)
#
# e <- marpoly_errors(fit)
# marpoly_diagnostics(e, powers = c(1, 2, 3, 4), lag.max = 20)
