
# Fitting of the univariate FPCs ------------------------------------------


# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


# Always
library(survival)
library(bamlss)
library(MFPCA)
library(sparseFLMM)
library(tidyverse)
library(parallel)
library(Rcpp)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/MJM_opt.R")
source("R/opt_updating_cpp.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing_cpp.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/fpca.R")
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")



# Simulate the data -------------------------------------------------------

# Prepare the model using the true FPCs
set.seed(1808)
mfpca <- create_true_MFPCA(M = 3, nmarker = 6, argvals = seq(0, 1, by = 0.01),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1,
                           evals = c(0.8, 0.5, 0.3))

# Simulate the data
simdat <- simMultiJM(nsub = 150, times = seq(0, 1, by = 0.01), 
                     probmiss = 0.75, maxfac = 2,
                     nmark = 6, long_assoc = "FPC", M = 6, 
                     FPC_bases = mfpca$functions, 
                     FPC_evals = mfpca$values,
                     mfpc_args = list(type = "split", eFunType = "PolyHigh",
                                      ignoreDeg = 1, eValType = "linear",
                                      eValScale = 1),
                     ncovar = 2,
                     lambda = function(t, x) {
                       1.6 * t^(0.6)
                     },
                     gamma = function(x) {
                       - 2 + 0.48*x[, 3]
                     },
                     alpha = list(function(t, x) {
                       0.64 + 0*t
                     }, function(t, x) {
                       -0.64 + 0*t
                     }, function(t, x) {
                       0.78 + 0*t
                     }, function(t, x) {
                       -0.78 + 0*t
                     }, function(t, x) {
                       0.93 + 0*t
                     }, function(t, x) {
                       -0.93 + 0*t
                     }),
                     mu = list(function(t, x, r){
                       2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }, function(t, x, r){
                       2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }, function(t, x, r){
                       1.5 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }, function(t, x, r){
                       1.5 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }, function(t, x, r){
                       1.1 + 0.3*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }, function(t, x, r){
                       1.1 + 0.23*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                     }),
                     sigma = function(t, x) {
                       log(0.2) + 0*t
                     }, 
                     tmax = NULL, seed = 1808, 
                     full = TRUE, file = NULL)


# Preprocessing: Estimate FPCs --------------------------------------------

# Prepare the model using estimated FPCs
set.seed(1808)
mfpca_own <- preproc_MFPCA(simdat$data, 
                          uni_mean = "y ~ 1 + obstime + x3 + obstime:x3")
mfpca_fdapace <- preproc_MFPCA(simdat$data, 
                          uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                          method = "FPCA")
mfpca_fpca.sc <- preproc_MFPCA(simdat$data, 
                               uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                               method = "fpca.sc")
mfpca_PACE <- preproc_MFPCA(simdat$data, 
                               uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                               method = "PACE")

# Plot the different FPC estimates
plot_fpca <- function(number = 1) {
  ref <- extractObs(mfpca_own$functions, number)
  ref1 <- flipFuns(ref, extractObs(mfpca_fdapace$functions, number))
  ref2 <- flipFuns(ref, extractObs(mfpca_fpca.sc$functions, number))
  ref3 <- flipFuns(ref, extractObs(mfpca_PACE$functions, number))
  dat <- data.frame(t = rep(unlist(argvals(ref)), 4),
                    y = unlist(
                          lapply(list(ref, ref1, ref2, ref3), function (x) {
                            lapply(x, function (mark) {
                              drop(mark@X)
                            })
                          })),
                    marker = factor(rep(rep(names(ref), 
                               each = length(argvals(ref[[1]])[[1]])), 4)),
                    type = factor(rep(c("own", "fdapace", "fpca.sc", "PACE"),
                                      each = length(unlist(argvals(ref))))))
  p <- ggplot(data = dat, aes(x = t, y = y, colour = type)) +
    geom_line() +
    facet_wrap(~marker, scales = "free") +
    ggtitle(paste("FPC", number))
  p
}

plot_fpca(1)
plot_fpca(2)
plot_fpca(3)
