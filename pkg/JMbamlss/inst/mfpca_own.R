# Covariance matrix of the random effects
auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
cross <- matrix(rep(0.03, 4), ncol = 2)
cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
                1, 0, 1, 0.75, 0.5, 0,
                0.75, 1, 0, 1, 0.75, 0.5,
                0.5, 0.75, 1, 0, 1, 0.75,
                0, 0.5, 0.75, 1, 0, 1,
                0, 0, 0.5, 0.75, 1, 0),
              ncol = 6)
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)), auto)

# MultiFunData object containing the basis functions
seq1 <- seq(0, 1, by = 0.01)
mfun1 <- multiFunData(
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                      byrow = TRUE, ncol = length(seq1))),
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                     byrow = TRUE, ncol = length(seq1))),
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                     byrow = TRUE, ncol = length(seq1))),
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                     byrow = TRUE, ncol = length(seq1))),
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                     byrow = TRUE, ncol = length(seq1))),
  funData(argvals = seq1,
          X = matrix(c(rep(1, length(seq1)), seq1),
                     byrow = TRUE, ncol = length(seq1)))
)

Basis_fun <- funData(argvals = seq1,
                     X = matrix(c(rep(1, length(seq1)), seq1),
                                byrow = TRUE, ncol = length(seq1)))
B_block <- MFPCA:::calcBasisIntegrals(Basis_fun@X, 1, Basis_fun@argvals)

Bchol <- Matrix::bdiag(rep(list(Matrix::chol(B_block)), 6))

vectors <- Re(eigen(Matrix::crossprod(Bchol) %*% cov)$vectors[, seq_len(12)])
