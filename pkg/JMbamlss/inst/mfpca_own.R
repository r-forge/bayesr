library(Matrix)
library(funData)
library(foreach)

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
Basis_fun <- funData(argvals = seq1,
                     X = matrix(c(rep(1, length(seq1)), seq1),
                                byrow = TRUE, ncol = length(seq1)))
B_block <- Matrix::bdiag(rep(list(
  MFPCA:::calcBasisIntegrals(Basis_fun@X, 1, Basis_fun@argvals)),
  6))

eig <- eigen(B_block %*% cov)
values <- Re(eig$values[seq_len(12)])
vectors <- Re(eig$vectors[, seq_len(12)])

npc <- c(2, 2, 2, 2, 2, 2)
npcCum <- cumsum(c(0, npc))

# normalization factors
normFactors <- 1/sqrt(diag(t(vectors) %*% cov %*% vectors))

tmpWeights <- cov %*% vectors
eFunctions <- multiFunData(foreach::foreach(j = seq_len(6)) %do% {
  MFPCA:::univExpansion(type = "given",
                scores = 1/sqrt(values) * normFactors * 
                  t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                argvals = list(seq1),
                functions = Basis_fun,
                params = NULL)
})

norm(eFunctions)


basis_funs <- rep(list(Basis_fun), 6)

#' Function to calculate the multivariate FPCs for a given covariance matrix
#' 
#' @param cov Covariance matrix of the basis functions coefficients.
#' @param basis_funs List with basis functions on each dimension. The basis
#'  functions are funData objects
MFPCA_cov <- function (cov, basis_funs) {
  
  # Information about the objects
  p <- length(basis_funs)
  npc <- sapply(basis_funs, nObs)
  npcCum <- cumsum(c(0, npc))
  argvals_list <- lapply(basis_funs, function(x) {
    argvals(x)[[1]]
  })
  
  # Construct matrix of basis function integrals
  B_block <- Matrix::bdiag(lapply(basis_funs, function (x){
    MFPCA:::calcBasisIntegrals(x@X, dimSupp = 1, x@argvals)
  }))
  
  # Eigendecomposition
  eig <- eigen(B_block %*% cov)
  values <- Re(eig$values)
  vectors <- Re(eig$vectors)
  
  # Components of estimation formula
  # before the sums
  normFactors <- 1/sqrt(diag(t(vectors) %*% cov %*% vectors))
  # after the sums
  blockWeights <- cov %*% vectors
  
  # Calculation of MFPCs
  eFunctions <- multiFunData(
    foreach:::'%do%'(foreach::foreach(j = seq_len(p)), {
      MFPCA:::univExpansion(type = "given",
                            scores = 1/sqrt(values) * normFactors * 
                              t(blockWeights[npcCum[j] + seq_len(npc[j]), , 
                                             drop = FALSE]),
                            argvals = argvals_list[j],
                            functions = basis_funs[[j]],
                            params = NULL)
  }))
  
  out <- list("values" = values,
              "functions" = eFunctions,
              "vectors" = vectors,
              "normFactors" = normFactors)
}

juhu <- MFPCA_cov(cov = cov, basis_funs = basis_funs)
