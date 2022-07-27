
# Set Up Small Simulation -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
location <- "server_linux"
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
library(tidyverse)
library(parallel)
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


# Setting for the simulation
start <- 100
stop <- 119
number_cores <- 2
setting <- "B"
dir.create(paste0("simulation/", setting), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Objects for all clusters ------------------------------------------------


# Prepare the model using the true FPCs
set.seed(1808)
mfpca <- create_true_MFPCA(M = 3, nmarker = 6, argvals = seq(0, 1, by = 0.01),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1,
                           evals = c(0.8, 0.5, 0.3))
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

simulate_models <- function (i){
  seed <- i
  
  # Simulate the data
  simdat <- simMultiJM(nsub = 150, times = seq(0, 1, by = 0.01), 
                       probmiss = 0.75, maxfac = 2,
                       nmark = 6, param_assoc = FALSE, M = 6, 
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
                       tmax = NULL, seed = seed, 
                       full = TRUE, file = NULL)
  # ggplot(data = simdat$data_short %>% slice_head(n = 150),
  #        aes(x = id, y = survtime, colour = as.factor(event)))+ geom_point()
  # ggplot(simdat$data, aes(x = obstime, y = y, color = id)) +
  #   geom_line() +
  #   facet_grid(~marker) +
  #   theme(legend.position = "none")
  
  
  # Prepare the model using estimated FPCs
  mfpca_es <- preproc_MFPCA(simdat$data, 
                            uni_mean = "y ~ 1 + obstime + x3 + obstime:x3")
  vals <- which(mfpca_es$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_es$values[vals])/sum(mfpca_es$values[vals]) > 0.95))
  mfpca_es_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_es) {
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
        paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
               "mfpca_es_list[[", x, "]]))")
      }), collapse = " + "))),
    sigma ~ -1 + marker,
    alpha ~ -1 + marker
  )
  
  # Estimate the model using the true FPCs
  b_tru <- try(bamlss(f_tru, family = mjm_bamlss, data = simdat$data, 
                      timevar = "obstime", maxit = 1500, n.iter = 23000,
                      burnin = 3000, thin = 20))
  
  # Estimate the model using the estimated FPCs
  b_est <- try(bamlss(f_est, family = mjm_bamlss, data = simdat_es, 
                      timevar = "obstime", maxit = 1500, n.iter = 23000,
                      burnin = 3000, thin = 20))
 
  # Output of simulation
  out <- list(simdat = simdat, simdat_es = simdat_es,
              b_tru = b_tru, b_est = b_est)
  try(save(out, file = paste0("simulation/", setting, "/sim", i, ".Rdata")))
  NULL
}



# Simulation --------------------------------------------------------------

# Actual parallel section
cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
cl <- makeForkCluster(number_cores)
sim <- parLapply(cl, start:stop, simulate_models)
stopCluster(cl)
Sys.time()

rm(list =ls())
quit("no")


 