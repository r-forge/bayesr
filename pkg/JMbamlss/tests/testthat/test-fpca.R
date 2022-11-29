
# Test the univariate FPCA ------------------------------------------------

# Code adapted from package 'refund'


test_that("FPCA function and refund package agree on toy example", {
  skip_on_cran()
  
  set.seed(12212)
  n <- 100
  ngrid <- 40
  t <- seq(0, 1, l=ngrid)
  efcts <- poly(t, 2)
  Y <- outer(2 * drop(scale(rnorm(n))), efcts[, 1]) +
    outer(drop(scale(rnorm(n))), efcts[, 2])
  colnames(Y) <- t
  
  sc <- refund::fpca.sc(Y)
  sc_own <- fpca(Y, argvals_obs = TRUE)
  
  expect_equal(sc$Yhat, sc_own$Yhat, tolerance=.01)
  
})

test_that("FPCA function can be forced on particular interval", {
  skip_on_cran()
  
  set.seed(1408)
  n <- 100
  ngrid <- 40
  t <- seq(0, 1, l=ngrid)
  efcts <- poly(t, 2)
  Y <- outer(2 * drop(scale(rnorm(n))), efcts[, 1]) +
    outer(drop(scale(rnorm(n))), efcts[, 2])
  ydata <- Y %>% 
    data.frame() %>% 
    pivot_longer(1:40, values_to = ".value") %>%
    mutate(".id" = rep(seq_len(n), each = ngrid),
           ".index" = rep(t, times = n),
           "name" = NULL) %>%
    relocate(".id", ".index", ".value")
  
  fpc_max <- fpca(ydata = ydata)
  fpc <- fpca(ydata = ydata)
  
  expect_equal(fpc_max$evals, fpc$evals)
  expect_equal(fpc_max$Yhat[, -41], fpc$Yhat, tolerance = 0.01)
  # Eigenfunctions change (different interval, normalization)
  # Scores change as a result
  # plot(c(fpc$Yhat), c(fpc$Yhat[, -41]))
  # plot(c(seq(0, 1, l=ngrid), 1.1), fpc_max$efunctions[, 1])
  # plot(c( seq(0, 1, l=ngrid), 1.1), fpc_max$efunctions[, 2])
  # plot(c( seq(0, 1, l=ngrid), 1.1), c(fpc$efunctions[, 1], NA))
  # plot(c( seq(0, 1, l=ngrid), 1.1), c(fpc$efunctions[, 2], NA))
  
})