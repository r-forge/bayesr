
# Set Up Small Simulation -------------------------------------------------



# Set up R session --------------------------------------------------------



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
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
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
source("R/mfpca_sim.R")
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")


# Setting for the simulation
start <- 100
stop <- 149
number_cores <- 2
setting <- "scen_I_130922"
Sys.time()
sessionInfo()

# Calculate the true multivariate FPCs

# Covariance matrix for the data generation
auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
cross <- matrix(rep(0.03, 4), ncol = 2)
cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
                1, 0, 1, 0.75, 0.5, 0,
                0.75, 1, 0, 1, 0.75, 0.5,
                0.5, 0.75, 1, 0, 1, 0.75,
                0, 0.5, 0.75, 1, 0, 1,
                0, 0, 0.5, 0.75, 1, 0),
              ncol = 6)
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)), 
                                         auto)

# Basis functions on each dimension
seq1 <- seq(0, 1, by = 0.01)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)


# Prepare the model using the true FPCs -----------------------------------

# Prepare objects for the model on different data sets
mfpca_tru <- MFPCA_cov(cov = cov, basis_funs = b_funs)

nfpc <- min(which(
  cumsum(mfpca_tru$values)/sum(mfpca_tru$values) > 0.95))
mfpca_tru_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_tru_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)



# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  load(paste0("simulation/", setting, "/data/d", i, ".Rdata"))

  try_obj <- try({
    
    # Estimate the model using tru FPCs
    d_rirs_tru <- attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
    t_tru <- system.time(
      b_est <- bamlss(f_tru, family = mjm_bamlss, data = d_rirs_tru, 
                      timevar = "obstime", maxit = 1500, n.iter = 5500,
                      burnin = 500, thin = 5)
    )
    attr(b_est, "comp_time") <- t_tru
    attr(b_est, "FPCs") <- mfpca_tru
    attr(b_est, "nfpc") <- nfpc
    save(b_est, file = paste0("simulation/", setting, "/bamlss_tru/b", i, 
                              ".Rdata"))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
