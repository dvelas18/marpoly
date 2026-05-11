test_that("stationarity check works", {
  st <- marpoly_stationary(causal = 0.7, noncausal = 0.2)
  expect_true(st$stationary)
  expect_error(marpoly_simulate(50, causal = 1.1, noncausal = NULL), "Stationarity")
})

test_that("simulation returns requested length", {
  sim <- marpoly_simulate(100, causal = 0.2, noncausal = 0.7, seed = 123)
  expect_equal(length(sim$x), 100)
})

test_that("periodogram dimensions are correct", {
  sim <- marpoly_simulate(100, causal = 0.2, noncausal = 0.7, seed = 123)
  p2 <- marpoly_periodogram(sim$x, order = 2, grid_size = 20)
  p3 <- marpoly_periodogram(sim$x, order = 3, grid_size = 10)
  p4 <- marpoly_periodogram(sim$x, order = 4, grid_size = 6)
  expect_equal(length(p2$values), length(p2$frequency))
  expect_equal(dim(p3$values), c(10, 10))
  expect_equal(dim(p4$values), c(6, 6, 6))
})

test_that("root allocations enumerate AR(2) possibilities", {
  allocs <- marpoly_root_allocations(c(0.9, -0.14))
  expect_equal(nrow(allocs), 4)
  expect_true(any(allocs$model == "AR(2,0)"))
  expect_true(any(allocs$model == "AR(0,2)"))
})
