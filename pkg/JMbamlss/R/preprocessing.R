
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
#' @param marker String giving the name of the longitudinal marker variable.
#' @param M Number of mFPCs to compute in the MFPCA. If not supplied, it
#'  defaults to the maximum number of computable mFPCs.
#' @param remove_obs Minimal number of observations per individual and marker to
#'  be included in the FPC estimation. Defaults to NULL (all observations). Not
#'  removing observations can lead to problems if the univariate variance 
#'  estimate is negative and has to be truncated, then the scores for IDs with
#'  few observations cannot be estimated.
#' @param method Which package to use for the univariate FPCA. Either function
#'  adapted function 'fpca', 'FPCA' from package \code{fdapace}, 'fpca.sc' from
#'  package \code{refund}, or function 'PACE' from package \code{MFPCA}.
#' @param var_names Vector of variable names corresponding to ID and MARKER. 
#'  Needed to remove observations.
#' @param nbasis Number of B-spline basis functions for mean estimate and
#'  bivariate smoothing of covariance surface for methods fpca.sc, PACE.
#' @param npc Number of univariate principal components to use in fpca.sc, PACE.
#' @param fve_uni Fraction of univariate variance explained for method FPCA.
#' @param pve_uni Proportion of univariate variance explained for method 
#'  fpca.sc,PACE.
#' @param save_uniFPCA TRUE to attach list of univariate FPCAs as attribute to
#'  output. Defaults to FALSE.
preproc_MFPCA <- function (data, uni_mean = "y ~ s(obstime) + s(x2)", 
                           time = "obstime", id = "id", marker = "marker",
                           M = NULL, remove_obs = NULL, 
                           method = c("fpca", "fpca.sc", "FPCA", "PACE"), 
                           nbasis = 10, npc = NULL,
                           fve_uni = 0.99, pve_uni = 0.99,
                           save_uniFPCA = FALSE) {
  require(bamlss)
  require(MFPCA)
  method <- match.arg(method)
  
  if (!is.null(remove_obs)) {
    few_obs <- apply(table(data[, id], data[, marker]), 1,
                     function (x) any(x < remove_obs))
    data <- droplevels(data[data[, id] %in% paste(which(!few_obs)), ])
    
  }
  
  marker_dat <- split(data, data$marker)
  
  # Check whether enough observations are available on each marker to be able
  # to estimate the same full interval on all markers
  maxtime <- sapply(marker_dat, function(x) max(x[, time]))
  if (length(unique(maxtime)) > 1) {
    if (method != "fpca") {
      stop("Estimation of MFPCA for different univariate time intervals ",
           "has only been tested for method 'fpca'. Note that the design ",
           "matrices constructed in the joint model will not be correct.")
    } else {
      warning("Estimation of univariate FPCAs has to extrapolate for ",
              paste(paste0("marker", which(maxtime != max(maxtime))),
                    sep = ", "),
              " from ", 
              paste(maxtime[which(maxtime != max(maxtime))], sep = ","), 
              " to ", max(maxtime), ". Please check if appropriate.")
    }
  }
  argvals_pred <- seq(0, max(maxtime), length.out = 101)
  
  uni_mean <- as.formula(uni_mean)
  
  marker_dat <- lapply(marker_dat, function (mark) {
    mark$res <- bam(formula = uni_mean, data = mark)$residuals
    mark
  })
  
  if (method == "fpca.sc" | method == "fpca") {
    require(refund)
    
    # Construct objects for fpca.sc function
    lY <- lapply(marker_dat, function (mark) {
      data.frame(".id" = mark[, id],
                 ".index" = mark[, time],
                 ".value" = mark$res)
    })
    
    # FPCA for each marker
    if (method == "fpca.sc") {
      FPCA <- lapply(lY, function(y) {
        refund::fpca.sc(ydata = y, pve = pve_uni, nbasis = nbasis, npc = npc,
                        var = TRUE)
      })
    } else {
      FPCA <- lapply(lY, function(y) {
        fpca(ydata = y, pve = pve_uni, nbasis = nbasis, npc = npc, var = TRUE,
             argvals_pred = argvals_pred)
      })
    }
    
    # Construct multivariate FunData and estimated FPCs
    mFData <- multiFunData(lapply(FPCA, function (mark) {
      funData(argvals = mark$argvals, X = mark$Yhat)
    }))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = funData(argvals = mark$argvals,
                                               X = t(mark$efunctions)))
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "npc"))
    }
    
    
  } else if (method == "FPCA") {
    require(fdapace)
    
    # Construct objects for FPCA function
    ly <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      split(mark$res, mark[, id])
    })
    lt <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      split(mark[, time], mark[, id])
    })
    
    # FPCA for each marker
    FPCA <- mapply(fdapace::FPCA, Ly = ly, Lt = lt, SIMPLIFY = FALSE,
                   MoreArgs = list(optns = list(FVEthreshold = fve_uni,
                                                nRegGrid = 101)))
    
    # Construct multivariate FunData and estimated FPCs
    mFData <- multiFunData(lapply(FPCA, function (mark) {
      funData(argvals = mark$workGrid, X = fitted(mark))
    }))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = funData(argvals = mark$workGrid,
                                               X = t(mark$phi)))
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "selectK"))
    }
    
  } else if (method == "PACE") {
    
    # Construct irregular FunData
    m_irregFunData <- lapply(marker_dat, function (mark) {
      mark <- mark[order(mark[, time]), ]
      irregFunData(argvals = split(mark[, time], mark[, id]), 
                   X = split(mark$res, mark[, id]))
    })
    
    # Remove observations with too few scalar observations
    if (!is.null(npc)) {
      rem <- lapply(m_irregFunData, function (mark) {
        which(lapply(mark@argvals, length) < npc)
      })
      rem <- Reduce(union, rem)
      take <- if (length(rem) == 0) {
        seq_len(nObs(m_irregFunData[[1]]))
      } else {
        seq_len(nObs(m_irregFunData[[1]]))[-rem]
      }
      m_irregFunData <- lapply(m_irregFunData, function (mark) {
        extractObs(mark, obs = take)
      })
    }
    
    # Use PACE function for each marker
    FPCA <- lapply(m_irregFunData, function(mark) {
      PACE(mark, npc = npc, pve = pve_uni, nbasis = nbasis)
    })
    
    # Construct multivariate FunData and estimated FPCs
    mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
    uniExpansions <- lapply(FPCA, function (mark) {
      list(type = "given", functions = mark$functions)
    })
    if (is.null(M)) {
      M <- sum(sapply(FPCA, "[[", "npc"))
    }
  }
  

  MFPCA <- MFPCA(mFData = mFData, M = M, uniExpansions = uniExpansions)
  attr(MFPCA, "sigma2") <- lapply(FPCA, "[[", "sigma2")
  if (save_uniFPCA) {
    attr(MFPCA, "uniFPCA") <- FPCA
  }
  MFPCA
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



# Attach Weighted Functional Principal Components to the Data -----------


#' Attach Weighted Functional Principal Components to the Data
#' 
#' @param mfpc MFPCA object from which to extract the weighted FPCS.
#' @param data Data set to which the weighted FPCS are to be attached.
#' @param n Number of FPCs to attach. Defaults to NULL which corresponds to all
#'  FPCs in mfpc.
#' @param obstime Name of the time variable in data set at which points to 
#'  evaluate.
#' @param marker Name of the marker variable in the data set which separates the
#'  data.
attach_wfpc <- function(mfpca, data, n = NULL, obstime = "obstime", marker = "marker",
                        eval_weight = FALSE){
  
  # Is the data sorted by marker
  if (!all(order(data[[marker]]) == seq_len(nrow(data)))) message("ORDER!")
  
  wfpc <- NULL
  
  splitdat <- split(data[[obstime]], data[[marker]])
  
  if (!is.null(n)) {
    mfpca <- list(values = mfpca$values[seq_len(n)],
                  functions = extractObs(mfpca$functions, obs = seq_len(n)))
  }
  
  # eval_mpfc evaluates on all markers, so choose only the current one
  for (mark in seq_along(levels(data[[marker]]))) {
    mobs <- length(splitdat[[mark]])
    tot_wfpc <- eval_mfpc(mfpca = mfpca, timepoints = splitdat[[mark]],
                          eval_weight = eval_weight)
    wfpc <- rbind(wfpc, tot_wfpc[(mark-1)*mobs + seq_len(mobs), , drop = FALSE])
  }
  
  data <- cbind(data, wfpc)
  data
  
}
