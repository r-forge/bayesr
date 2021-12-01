
# Preprocessing Steps for Model Fitting -----------------------------------



# MFPCA object from data --------------------------------------------------

#' Preprocessing step to create MFPCA object
#' 
#' This function takes the data und uses the residuals of marker-specific
#' additive models to estimate the covariance structure for a MFPCA
#' 
#' @param data Data.frame such as returned by function simMultiJM.
#' @param uni_mean String to crate a formula for the univariate addtive models.
#' @param time String giving the name of the longitudinal time variable.
#' @param id String giving the name of the identifier.
#' @param M Number of mFPCs to compute in the MFPCA.
preproc_MFPCA <- function (data, uni_mean = "y ~ s(obstime) + s(x2)", 
                           time = "obstime", id = "id", M = 2) {
  require(bamlss)
  require(MFPCA)
  
  marker_dat <- split(data, data$marker)
  
  uni_mean <- as.formula(uni_mean)
  
  marker_dat <- lapply(marker_dat, function (mark) {
    mark$res <- bam(formula = uni_mean, data = mark)$residuals
    mark
  })
  m_irregFunData <- lapply(marker_dat, function (mark) {
    mark <- mark[order(mark[, time]), ]
    irregFunData(argvals = split(mark[, time], mark[, id]), 
                 X = split(mark$res, mark[, id]))
  })
  FPCA <- lapply(m_irregFunData, function(mark) {
    PACE(mark)
  })
  mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
  uniExpansions <- lapply(FPCA, function (mark) {
    list(type = "given", functions = mark$functions)
  })
  MFPCA <- MFPCA(mFData = mFData, M = M, uniExpansions = uniExpansions)
}



# Create true MFPC basis --------------------------------------------------

mfpc_args = list(type = "split", eFunType = "Poly",
                 ignoreDeg = NULL, eValType = "linear",
                 eValScale = 1)

create_true_MFPCA <- function (M, nmarker, argvals = seq(0, 120, 1), 
                               type = "split", eFunType = "Poly",
                               ignoreDeg = NULL, eValType = "linear",
                               eValScale = 1) {
  
  evals <- funData::eVal(M = M, type = eValType)
  evals <- eValScale * evals
  
  mfpc_seed <- switch(type, 
                      "split" = sample(c(-1, 1), nmarker, 0.5),
                      "weight" = stats::runif(nmark, 0.2, 0.8))
  
  bases <- switch(type,
                  split = simMultiSplit(argvals = argvals, M = M,
                                        eFunType = eFunType,
                                        ignoreDeg = ignoreDeg,
                                        eValType = eValType,
                                        s = mfpc_seed),
                  weighted = simMultiWeight(argvals = argvals, M = M,
                                            eFunType = eFunType,
                                            ignoreDeg = ignoreDeg,
                                            eValType = eValType,
                                            alpha = mfpc_seed),
                  stop(paste0("Choose either 'split' or 'weighted' for the sim",
                              "ulation of multivariate functional data.")))
  
}
