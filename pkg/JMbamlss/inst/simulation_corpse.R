
# Set Up Small Simulation -------------------------------------------------


# Set up R session --------------------------------------------------------



# Specify location
location <- "laptop"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") "H:/volkmana.hub/R4_linux"
            else "H:/R4_windows")
  setwd(if (location == "server_linux") "H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


# Always
library(survival)
library(bamlss)
library(MFPCA)
library(tidyverse)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/compile.R")
compile_alex(location)




# Objects for all clusters ------------------------------------------------


# Prepare the model using the true FPCs
mfpca <- create_true_MFPCA(M = 3, nmarker = 2, argvals = seq(0, 1, by = 0.01),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(0.8, 0.5, 0.3))
mfpca_list <- lapply(1:3, function (i, m = mfpca) {
  list(functions = extractObs(m$functions, i),
       values = m$values[i])
})
f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + 
    s(id, fpc.1, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[1]])) + 
    s(id, fpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[2]])) + 
    s(id, fpc.3, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[3]])),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)



# Simulation function -----------------------------------------------------

simulate_models <- function (i, start){
  seed <- (start + i)
  
  # Simulate the data
  simdat <- simMultiJM(nsub = 150, times = seq(0, 1, by = 0.01), 
                       probmiss = 0.75, maxfac = 2,
                       nmark = 2, param_assoc = FALSE, M = 3, 
                       FPC_bases = NULL, 
                       FPC_evals = c(0.8, 0.5, 0.3),
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
                       }),
                       mu = list(function(t, x, r){
                         2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                       }, function(t, x, r){
                         2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                       }),
                       sigma = function(t, x) {
                         log(0.2) + 0*t
                       }, 
                       tmax = NULL, seed = seed, 
                       full = TRUE, file = NULL)
  
  # Prepare the model using estimated FPCs
  mfpca_es <- preproc_MFPCA(simdat$data, 
                            uni_mean = "y ~ 1 + obstime + x3 + obstime:x3")
  nfpc <- min(which(cumsum(mfpca_es$values)/sum(mfpca_es$values) > 0.95))
  mfpca_es_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_es) {
    list(functions = extractObs(mfpca$functions, i),
         values = mfpca$values[i])
  })
  simdat_es <- subset(simdat$data, 
                      select = -grep("fpc\\.", colnames(simdat$data)))
  simdat_es <- attach_wfpc(mfpca_es, simdat_es, n = nfpc)
  f_est <- list(
    Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10, bs = "ps"),
    gamma ~ 1 + x3,
    as.formula(paste0(
      "mu ~ -1 + marker + obstime:marker + x3:marker + ",
      paste0(lapply(seq_len(nfpc), function(x) {
        paste0("s(id, fpc.", x, ", bs = unc_pcre, xt = list('mfpc' =",
               " mfpca_list[[", x, "]]))")
      }), collapse = " + "))),
    sigma ~ -1 + marker,
    alpha ~ -1 + marker
  )
  
}

}
 