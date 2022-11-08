test_that("survival integral for survival predictors works", {
  
  # Create some toy data
  pre_fac <- c(0.5, 1, 1.5)
  pre_vec <- matrix(c(2, 8), ncol = 2, nrow = 3, byrow = TRUE)
  int_vec_lambda <- matrix(c(0.2, 0.8, 0.4, 0.4, 0.1, 0.9,
                             0.1, 0.7, 0.3, 0.3, 0.1, 0.8,
                             0.3, 0.9, 0.5, 0.5, 0.2, 1), ncol = 2,
                           byrow = TRUE)
  omega <- rep(c(0.1, 0.5, 1), times = 3)
  weights <- statmod::gauss.quad(3)$weights
  survtime <- c(0.2, 0.7, 1)
  
  # Manual calculation of the integrals
  man_lambda <- {
    out <- list("score_int" = vector(length = 2), 
                "hess_int" = vector(length = 4))
    for (i in 1:3) {
      out$score_int <- out$score_int + pre_fac[i] * survtime[i]/2 *
        drop(weights %*% 
               (omega[(i-1)*3 + 1:3] * int_vec_lambda[(i-1)*3 + 1:3, ]))
      out$hess_int <- out$hess_int + pre_fac[i] * survtime[i]/2 * 
        drop(weights %*% 
                  t(sapply(seq_len(3), function (j) {
                    omega[(i-1)*3 + j] * 
                      c(tcrossprod(int_vec_lambda[(i-1)*3 + j, ]))
                  })))
    }
    out
  }
  
  man_gamma <- {
    out <- list("score_int" = vector(length = 2), 
                "hess_int" = vector(length = 4))
    for (i in 1:3) {
      out$score_int <- out$score_int + pre_fac[i] * pre_vec[i, ] * 
        survtime[i]/2 * drop(weights %*% omega[(i-1)*3 + 1:3])
      out$hess_int <- out$hess_int + pre_fac[i] * c(tcrossprod(pre_vec[i, ])) * 
        survtime[i]/2 * drop(weights %*% omega[(i-1)*3 + 1:3])
    }
    out
  }
  
  # C implementation of survival integral
  s_lambda <- survint_C(pred = "lambda", pre_fac, pre_vec = NULL, omega, 
                        int_fac = NULL, int_vec_lambda, weights, survtime)
  s_gamma <- survint_C(pred = "gamma", pre_fac, pre_vec, omega, int_fac = NULL,
                       int_vec = NULL, weights, survtime)
  
  # See if they align
  expect_equal(s_lambda, man_lambda)
  expect_equal(s_gamma, man_gamma)
})

test_that("Survival integral for longitudinal predictors work", {
  
  # Create some toy data
  # int_vec_long1: A design matrix as constructed by a marker-interaction
  # int_vec_long2: A design matrix with same coefficients over the markers
  # int_vec_pcre: A design matrix as constructed by a PCRE-term
  pre_fac <- c(0.5, 1, 1.5)
  int_vec_onedim <- matrix(c(0.2, 0.8, 0.4, 0.4, 0.1, 0.9,
                             0.1, 0.7, 0.3, 0.3, 0.1, 0.8,
                             0.3, 0.9, 0.5, 0.5, 0.2, 1), ncol = 2,
                           byrow = TRUE)
  int_vec_onepcre <- matrix(c(0.2, 0.4, 0.6, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0.1, 0.5, 0.7, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0.2, 0.4, 0.6), ncol = 3)
  int_vec_long1 <- rbind(int_vec_onedim, matrix(0, nrow = 9, ncol = 2))
  int_vec_long2 <- rbind(int_vec_onedim, 2*int_vec_onedim)
  int_vec_pcre <- rbind(int_vec_onepcre, 2*int_vec_onepcre)
  int_fac <- rep(c(0.5, 1), each = 3*3)
  omega <- rep(c(0.1, 0.5, 1), times = 3)
  weights <- statmod::gauss.quad(3)$weights
  survtime <- c(0.2, 0.7, 1)
  
  # Manual calculation of the integrals
  man_long1 <- {
    out <- list("score_int" = vector(length = 2), 
                "hess_int" = vector(length = 4))
    for (i in 1:3) {
      out$score_int <- out$score_int + pre_fac[i] * survtime[i]/2 *
        drop(weights %*% 
               (omega[(i-1)*3 + 1:3] * int_fac[(i-1)*3 + 1:3] * 
                  int_vec_long1[(i-1)*3 + 1:3, ] +
                  omega[(i-1)*3 + 1:3] * int_fac[(i-1)*3 + 1:3 + 9] *
                  int_vec_long1[(i-1)*3 + 1:3 + 9, ]))
      out$hess_int <- out$hess_int + pre_fac[i] * survtime[i]/2 * 
        (drop(weights %*% (
                t(sapply(seq_len(3), function (j) {
                  omega[(i-1)*3 + j] * int_fac[(i-1)*3 + j]^2 *
                    c(tcrossprod(int_vec_long1[(i-1)*3 + j, ]))
                  })) +
                  t(sapply(seq_len(3), function (j) {
                    omega[(i-1)*3 + j] * int_fac[(i-1)*3 + j + 9]^2 *
                      c(tcrossprod(int_vec_long1[(i-1)*3 + j + 9, ]))
                  })))))
    }
    out
  }
  
  man_long2 <- {
    out <- list("score_int" = vector(length = 2), 
                "hess_int" = vector(length = 4))
    for (i in 1:3) {
      out$score_int <- out$score_int + pre_fac[i] * survtime[i]/2 *
        drop(weights %*% 
               (omega[(i-1)*3 + 1:3] * (int_fac[(i-1)*3 + 1:3] * 
                  int_vec_long2[(i-1)*3 + 1:3, ] + 
                    int_fac[(i-1)*3 + 1:3 + 9] *
                  int_vec_long2[(i-1)*3 + 1:3 + 9, ])))
      out$hess_int <- out$hess_int + pre_fac[i] * survtime[i]/2 * 
        (drop(weights %*% (omega[(i-1)*3 + 1:3] *
          t(sapply(seq_len(3), function (j) {
            sum_vec <- int_fac[(i-1)*3 + j] * int_vec_long2[(i-1)*3 + j, ] +
              t(int_fac[(i-1)*3 + j + 9] * int_vec_long2[(i-1)*3 + j + 9, ])
            crossprod(sum_vec)
          })))))
    }
    out
  }
  # Use the same code as for man_long2 but with the by-marker design matrix
  # int_vec_long1: This gives the same result as man_long1
  
  man_pcre <- {
    out <- list("score_int" = vector(length = 3), 
                "hess_int" = vector(length = 9))
    for (i in 1:3) {
      out$score_int <- out$score_int + pre_fac[i] * survtime[i]/2 *
        drop(weights %*% 
               (omega[(i-1)*3 + 1:3] * (int_fac[(i-1)*3 + 1:3] * 
                                          int_vec_pcre[(i-1)*3 + 1:3, ] + 
                                          int_fac[(i-1)*3 + 1:3 + 9] *
                                          int_vec_pcre[(i-1)*3 + 1:3 + 9, ])))
      out$hess_int <- out$hess_int + pre_fac[i] * survtime[i]/2 * 
        (drop(weights %*% (omega[(i-1)*3 + 1:3] *
                             t(sapply(seq_len(3), function (j) {
                               sum_vec <- int_fac[(i-1)*3 + j] * 
                                 int_vec_pcre[(i-1)*3 + j, ] +
                                 t(int_fac[(i-1)*3 + j + 9] * 
                                     int_vec_pcre[(i-1)*3 + j + 9, ])
                               crossprod(sum_vec)
                             })))))
    }
    out
  }
  
  
  # Use C-implementation for calculation
  s_long1 <- survint_C(pred = "long", pre_fac, pre_vec = NULL, omega, int_fac, 
                       int_vec_long1, weights, survtime)
  s_long2 <- survint_C(pred = "long", pre_fac, pre_vec = NULL, omega, int_fac, 
                       int_vec_long2, weights, survtime)
  s_pcre <- survint_C(pred = "fpc_re", pre_fac, pre_vec = NULL, omega, int_fac,
                      int_vec_pcre, weights, survtime)
  s_pcre$hess_int <- c(diag(s_pcre$hess_int))
  s_pcre1 <- survint_C(pred = "long", pre_fac, pre_vec = NULL, omega, int_fac,
                       int_vec_pcre, weights, survtime)
  
  # Check whether manual and C implementation give the same results.
  expect_equal(s_long1, man_long1)
  expect_equal(s_long2, man_long2)
  expect_equal(s_pcre, man_pcre)
  expect_equal(s_pcre1, man_pcre)
  
})
