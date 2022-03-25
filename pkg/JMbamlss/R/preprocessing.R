
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
#' @param npc Number of univariate principal components to use in PACE.
preproc_MFPCA <- function (data, uni_mean = "y ~ s(obstime) + s(x2)", 
                           time = "obstime", id = "id", M = 2, npc = NULL) {
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
  if (!is.null(npc)) {
    rem <- lapply(m_irregFunData, function (mark) {
      which(lapply(mark@argvals, length) < npc)
    })
    rem <- Reduce(union, rem)
    take <- seq_len(nObs(m_irregFunData[[1]]))[-rem]
    m_irregFunData <- lapply(m_irregFunData, function (mark) {
      extractObs(mark, obs = take)
    })
  }
  FPCA <- lapply(m_irregFunData, function(mark) {
    PACE(mark, npc = npc)
  })
  mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
  uniExpansions <- lapply(FPCA, function (mark) {
    list(type = "given", functions = mark$functions)
  })
  MFPCA <- MFPCA(mFData = mFData, M = M, uniExpansions = uniExpansions)
}



# Create true MFPC basis --------------------------------------------------


create_true_MFPCA <- function (M, nmarker, argvals = seq(0, 120, 1), 
                               type = "split", eFunType = "Poly",
                               ignoreDeg = NULL, eValType = "linear",
                               eValScale = 1, evals = NULL) {
  
  if (is.null(evals)) {
    evals <- funData::eVal(M = M, type = eValType)
    evals <- eValScale * evals
  }
  
  mfpc_seed <- switch(type, 
                      "split" = sample(c(-1, 1), nmarker, 0.5),
                      "weight" = stats::runif(nmark, 0.2, 0.8))
  
  mean <- multiFunData(
    rep(list(funData(argvals, matrix(0, nrow = M, ncol = length(argvals)))),
             nmarker))
  argvals <- rep(list(argvals), nmarker)
  
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
  
  if (M == 1) {
    bases <- multiFunData(lapply(bases, "/", sqrt(norm(bases))))
  }
  
  mfpca <- list(values = evals,
                functions = bases,
                meanFunction = mean)
  return(mfpca)
}


#------------------------------------------------------------------------------#
# Attach Weighted Functional Principal Components to the Data
#------------------------------------------------------------------------------#
#' Attach Weighted Functional Principal Components to the Data
#' 
#' @param mfpc MFPCA object from which to extract the weighted FPCS.
#' @param data Data set to which the weighted FPCS are to be attached.
attach_wfpc <- function(mfpca, data, obstime = "obstime", marker = "marker"){
  
  # Is the data sorted by marker
  if (all(order(data[[marker]]) != seq_len(nrow(data)))) message("ORDER!")
  
  wfpc <- NULL
  
  splitdat <- split(data[[obstime]], data[[marker]])
  
  # eval_mpfc evaluates on all markers, so choose only the current one
  for (mark in seq_along(levels(data[[marker]]))) {
    mobs <- length(splitdat[[mark]])
    tot_wfpc <- eval_mfpc(mfpca = mfpca, timepoints = splitdat[[mark]])
    wfpc <- rbind(wfpc, tot_wfpc[(mark-1)*mobs + seq_len(mobs), ])
  }
  
  data <- cbind(data, wfpc)
  data
  
}
