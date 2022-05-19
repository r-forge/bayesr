# 
# set.seed(1)
# 
# ### simulate data (one-dimensional domains)
# sim <-  simMultiFunData(type = "split", 
#                         argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#                         M = 5, eFunType = "Poly", eValType = "linear", 
#                         N = 100)
# 
# # MFPCA based on univariate FPCA
# uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#                                                         list(type = "uFPCA")))
# 
# fu_onen <- funData(argvals = seq(0, 1, length = 10),
#                    X = matrix(sin(seq(0, 3, length = 10)), nrow = 1))
# fu_muln <- funData(argvals = seq(0, 1, length = 10),
#                    X = matrix(c(sin(seq(0, 3, length = 10)),
#                                 cos(seq(0, 3, length = 10))),
#                               ncol = 10, nrow = 2, byrow = TRUE))

eval_fundata <- function(funData, evalpoints) {
  # Extract info
  args <- unlist(argvals(funData))
  arglength <- length(args)
  evalength <- length(evalpoints)
  indx <- seq_len(arglength + evalength)
  argfull <- c(args, evalpoints)
  order <- order(argfull)
  n <- nrow(X(funData))
  # Blow up matrix with NA values and order it
  na_mat <- matrix(NA, ncol = evalength, nrow = n)
  xobs <- cbind(X(funData), na_mat)[, order]
  # Different cases as cbind returns vector when n == 1
  if (n > 1) {
    eval <- zoo::na.approx(t(xobs), x = argfull[order], na.rm = FALSE)
    matrix(eval[match(arglength + seq_len(evalength), order), ], ncol = n,
           nrow = evalength)
  } else {
    eval <- zoo::na.approx(xobs, x = argfull[order], na.rm = FALSE)
    matrix(eval[match(arglength + seq_len(evalength), order)], ncol = n, 
           nrow = evalength)
  }
}

# eval_fundata(fu_onen, 0.39)
# eval_fundata(fu_muln, 0.39)
# eval_fundata(fu_onen, c(0.39, 0.7))
# eval_fundata(fu_muln, c(0.39, 0.7))


eval_mfpc <- function(mfpca, timepoints, marker = NULL, eval_weight = FALSE) {
  eigenfct <- mfpca$functions
  eigenval <- mfpca$values
  K <- length(eigenfct)
  
  if (is.null(marker)) {
    
    # predict the same time points for all markers
    prm <- do.call(rbind, 
                   lapply(eigenfct, eval_fundata, evalpoints = timepoints))
    
    rownames(prm) <- paste0("m", rep(seq_len(K), each = length(timepoints)))
  
  } else {
    
    # use different time points for markers
    if (length(timepoints) != length(marker)) {
      stop("Marker-specific prediction of the FPCs only possible ",
           "when the length of timepoints and marker values are equal.")
    }
    marker_time <- split(timepoints, marker)
    prm <- do.call(rbind, mapply(eval_fundata, funData = eigenfct, 
                                 evalpoints = marker_time, SIMPLIFY = FALSE))
    rownames(prm) <- unlist(lapply(seq_along(marker_time), function (x) {
      rep(paste0("m", x), length(marker_time[[x]]))
    }))
    
    # make sure that the ordering is correct
    ordering <- unlist(split(seq_along(timepoints), marker))
    prm <- prm[order(ordering), , drop = FALSE]
  }

  if (eval_weight) {
    colnames(prm) <- paste0("wfpc.", seq_len(ncol(prm)))
    # weight with eigenvalues
    t(t(prm) * sqrt(eigenval))
  } else {
    colnames(prm) <- paste0("fpc.", seq_len(ncol(prm)))
    # no weighting
    t(t(prm))
  }

}

# eval_mfpc(uFPCA, 0.39)
# eval_mfpc(uFPCA, c(0.39, 0.41))
# eval_mfpc(uFPCA, c(0.39, 0.7))
# eval_mfpc(uFPCA, c(0.39, 0.7, 0.39, 0.7), factor(c("m1", "m1", "m2", "m2")))
# eval_mfpc(uFPCA, c(0.39, 0.41, 0.39, 0.41), factor(c("m1", "m2", "m2", "m1")))
