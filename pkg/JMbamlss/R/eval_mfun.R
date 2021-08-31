
# set.seed(1)
# 
# ### simulate data (one-dimensional domains)
# sim <-  simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#                         M = 5, eFunType = "Poly", eValType = "linear", N = 100)
# 
# # MFPCA based on univariate FPCA
# uFPCA <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#                                                         list(type = "uFPCA")))


#' Evaluate 
#' @param MFPCA
#' @param timepoints
eval_mfpc <- function (MFPCA, timepoints) {
  
  # Which observations have to be interpolated
  # Keep track of the order to order the X values appropriately
  grid <- unlist(argvals(MFPCA$functions), recursive = FALSE)
  arg_vals <- mapply(function(tpoints, gpoints){
    unique(c(gpoints, tpoints))
  }, tpoints = timepoints, gpoints = grid, SIMPLIFY = FALSE)
  arg_order <- lapply(arg_vals, order)
  
  # Attach NAs to the observed values to later interpolate them
  obs_x <- lapply(MFPCA$functions, X)
  obs_x <- mapply(function(dim_x, dim_a){
    add_na <- cbind(dim_x, matrix(NA, nrow = nrow(dim_x),
                                  ncol = length(dim_a) - ncol(dim_x)))
    add_na
  }, dim_x = obs_x, dim_a = arg_vals, SIMPLIFY = FALSE)
  
  # Reorder the argvals and xvalues
  arg_vals <- mapply(function(tpoints, torder){
    tpoints[torder]
  }, tpoints = arg_vals, torder = arg_order, SIMPLIFY = FALSE)
  obs_x <- mapply(function(xpoints, xorder){
    xpoints[, xorder]
  }, xpoints = obs_x, xorder = arg_order, SIMPLIFY = FALSE)
  
  # linear interpolation
  na_evals <- mapply(function(xpoints, xvals) {
    zoo::na.approx(t(xpoints), x = xvals, na.rm = FALSE)
  }, xpoints = obs_x, xvals = arg_vals, SIMPLIFY = FALSE)

  # extract
  mapply(function (dim, tpoint) {
    rownames(dim) <- tpoint
    dim[paste(timepoints), ]
  }, dim = na_evals, tpoint = arg_vals, SIMPLIFY = FALSE)
}


# fu <- funData(argvals = seq(0, 1, length = 10), 
#               X = matrix(sin(seq(0, 3, length = 10)), nrow = 1))
# eval_fundata(fu, 0.39)

eval_fundata <- function(funData, evalpoints) {
  args <- unlist(argvals(funData))
  arglength <- length(args)
  evalength <- length(evalpoints)
  indx <- seq_len(arglength + evalength)
  argfull <- c(args, evalpoints)
  order <- order(argfull)
  xobs <- cbind(X(funData),
                matrix(NA, ncol = evalength, nrow = nrow(X(funData))))[, order]
  eval <- zoo::na.approx(t(xobs), x = argfull[order], na.rm = FALSE)
  eval[match(arglength + seq_len(evalength), order), ]
}
