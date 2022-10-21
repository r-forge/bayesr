#' Function to calculate the multivariate FPCA for a given covariance matrix and
#' univariate basis functions
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
