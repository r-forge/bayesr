test_that("Psi matrix multiplication works", {
  
  # Create some toy data
  d <- data.frame(X1 = c(0.5, 0.8, 1, 0, 0, 0.1, 0, 0, 0, 0.4, 0.5, 0),
                  X2 = c(0, 0, 0, 1, 0.4, 0, 0.3, 0.4, 0.1, 0, 0, 0.6),
                  id = factor(c(1, 1, 1, 2, 2, 1, 2, 2, 2, 1, 1, 2)))
  x <- list(X = cbind(d$X1, d$X2),
            cp_info = list(order(d$id),
                           as.integer(table(d$id))))
  y <- c(3, 2, 3, 4, 3, 1, 3, 2, 5, 2, 3, 1)
  mu <- c(1, 2, 3, 2, 3, 4, 1, 2, 3, 4, 2, 2)
  sigma <- rep(c(0.2, 0.4, 0.2), times = c(5, 4, 3))
  
  # C++ implementation
  H_cpp <- eigenMapMatMult(t(x$X * (1 / exp(sigma)^2)), x$X)
  
  # C implementation
  S_c <- psi_mat_crossprod(Psi = x, y = (y - mu) / exp(sigma)^2)
  H_c <- diag(psi_mat_crossprod(Psi = x, R = 1 / exp(sigma)^2)) 
  
  # R implementation
  S_r <- drop(crossprod(x$X, (y - mu) / exp(sigma)^2))
  H_r <- crossprod(x$X* (1 / exp(sigma)^2), x$X)
  
  expect_equal(H_c, H_r)
  expect_equal(H_cpp, H_r)
  expect_equal(S_c, S_r)
})
