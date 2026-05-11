# marpoly

`marpoly` provides R tools for **mixed causal-noncausal autoregressive models** and **polyspectral estimation**. The name combines MAR models with polyspectra: spectra, bispectra, and fourth-order spectra.

The package is designed for empirical users who want a direct workflow:

1. select the pseudo-causal AR order with second-order methods;
2. enumerate possible root allocations;
3. estimate AR(r,0), AR(0,s), and MAR(r,s) allocations with second-, third-, and fourth-order frequency-domain objectives;
4. rank allocations by the polyspectral objective;
5. extract model errors and check residual dependence.

## Installation

```r
# install.packages("remotes")
remotes::install_github("dvelas18/marpoly")
```

For local development:

```r
# install.packages("devtools")
devtools::load_all()
```

## Main functions

| Function | Purpose |
| --- | --- |
| `marpoly_simulate()` | Simulate stationary MAR(r,s) processes with normal, Student t, skew Student t, or user-defined innovations. |
| `marpoly_stationary()` | Check that causal and noncausal MAR polynomial roots are outside the unit circle. |
| `marpoly_periodogram()` | Compute order-2, order-3, and order-4 periodograms. |
| `marpoly_spectrum()` | Compute theoretical order-2, order-3, and order-4 MAR spectral densities. |
| `marpoly_plot()` | Plot a series with selected periodograms and fitted spectra. |
| `marpoly_ar_ols()` | Estimate the pseudo-causal AR(p) representation by OLS and report regression-style output. |
| `marpoly_select_p()` | Select the pseudo-causal AR order by AIC or BIC. |
| `marpoly_root_allocations()` | Enumerate all feasible root allocations for a pseudo-causal AR(p). |
| `marpoly_fit()` | Estimate and rank all selected root allocations using spectral orders 2, 3, and/or 4. |
| `marpoly_rank()` | Return a compact ranked table of root allocations and main competitors. |
| `marpoly_errors()` | Extract fitted errors for the selected allocation or for every allocation. |
| `marpoly_diagnostics()` | Check ACF and Ljung-Box diagnostics for residual powers. |

## Model convention

The package uses

```text
(1 - causal_1 L - ... - causal_r L^r)
(1 - noncausal_1 L^-1 - ... - noncausal_s L^-s) y_t = epsilon_t.
```

A stationary MAR(r,s) process requires both component polynomials to have roots outside the unit circle. `marpoly_simulate()` checks this restriction and stops with an error if it is violated.

## Quick start: simulate and estimate a MAR(1,1)

```r
library(marpoly)

sim <- marpoly_simulate(
  n = 600,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 20,
  gamma = 1.5,
  seed = 123
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
```

The summary reports the selected root allocation in a regression-like table with coefficient estimates, numerical-Hessian standard errors, test statistics, and p-values. The rank table reports the selected model and the main competing allocations without overwhelming the report with every diagnostic.

## Periodograms, spectra, and plots

```r
sim <- marpoly_simulate(
  n = 500,
  causal = 0.25,
  noncausal = 0.65,
  dist = "skewt",
  df = 8,
  gamma = 1.4,
  seed = 42
)

p2 <- marpoly_periodogram(sim$x, order = 2)
p3 <- marpoly_periodogram(sim$x, order = 3, grid_size = 32)
p4 <- marpoly_periodogram(sim$x, order = 4, grid_size = 12)

plot(p2)
plot(p3, component = "real")
plot(p4, component = "modulus", slice = 4)

s3 <- marpoly_spectrum(
  causal = 0.25,
  noncausal = 0.65,
  order = 3,
  kappa = 1,
  grid_size = 32
)
plot(s3, component = "real")

marpoly_plot(sim$x, orders = c(2, 3), grid_size = 32)
```

## Identification examples

### MAR(1,1) using orders 2 and 3 with skew Student t innovations

```r
sim <- marpoly_simulate(
  n = 800,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 20,
  gamma = 1.5,
  seed = 1
)
fit23 <- marpoly_fit(sim$x, p = 2, orders = c(2, 3), grid_size = 24)
summary(fit23)
```

### MAR(1,1) using orders 2 and 4 with Student t innovations

```r
sim <- marpoly_simulate(
  n = 800,
  causal = 0.20,
  noncausal = 0.70,
  dist = "t",
  df = 5,
  seed = 2
)
fit24 <- marpoly_fit(sim$x, p = 2, orders = c(2, 4), grid_size = 12)
summary(fit24)
```

### MAR(1,1) using orders 2, 3, and 4 with skew Student t innovations

```r
sim <- marpoly_simulate(
  n = 800,
  causal = 0.20,
  noncausal = 0.70,
  dist = "skewt",
  df = 5,
  gamma = 1.5,
  seed = 3
)
fit234 <- marpoly_fit(sim$x, p = 2, orders = c(2, 3, 4), grid_size = 12)
summary(fit234)
```

## Extracting errors for each allocation

```r
errors_all <- marpoly_errors(fit23, allocation = "all")

# ACF diagnostics for the selected allocation
selected_errors <- marpoly_errors(fit23)
diag_selected <- marpoly_diagnostics(selected_errors, powers = c(1, 2), lag.max = 20)
print(diag_selected)

# Compare allocations by Ljung-Box p-values for errors and squared errors
allocation_checks <- lapply(errors_all, marpoly_diagnostics, powers = c(1, 2), lag.max = 20)
```

## Additional examples included

The `examples/` folder contains scripts for:

- creating periodograms, spectral densities, and plots;
- identifying MAR(1,1) with orders 2 and 3 under skew Student t innovations;
- identifying MAR(1,1) with orders 2 and 4 under Student t innovations;
- identifying MAR(1,1) with orders 2, 3, and 4;
- estimation and root selection for MAR(1,1);
- estimation and root selection for MAR(1,2);
- extracting errors for each root allocation and checking the ACF of errors and squared errors;
- sensitivity to the frequency grid and order set;
- a real-data template.

## Practical notes

- Use `orders = c(2, 3)` when skewness is informative.
- Use `orders = c(2, 4)` when the innovation distribution is symmetric but heavy-tailed.
- Use `orders = c(2, 3, 4)` when both skewness and kurtosis may matter.
- Fourth-order estimation is more computationally demanding because it works on a three-dimensional frequency grid. Start with `grid_size = 10` or `grid_size = 12` and increase it for sensitivity checks.
- The selected root allocation is the one with the lowest objective value. The full allocation table is stored in `fit$results`, while `marpoly_rank(fit)` gives a compact user-facing ranking.

## References

Hecq, A. and Velásquez-Gaviria, D. (2025). **Spectral estimation for mixed causal-noncausal autoregressive models**. *Econometric Reviews*, 44(7), 939–962. DOI: `10.1080/07474938.2025.2465372`.

Velasco, C. and Lobato, I. N. (2018). **Frequency domain minimum distance inference for possibly noninvertible and noncausal ARMA models**. *The Annals of Statistics*, 46(2), 555–579. DOI: `10.1214/17-AOS1560`.
