
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
source("R/compile.R")
source("R/mfpca_sim.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")

# Load the data
load("simulation/scen_I_130922/data/d110.Rdata")


# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                 function (x) any(x < 2))
mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                             filter(id %in% paste(which(!few_obs))) %>% 
                             droplevels(), 
                           uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                           npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
d_rirs_est <- attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
f_long <- list(
  as.formula(paste0(
    "y ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker
)

# Model fit
t_est <- system.time(
  b_est <- bamlss(f_long, data = d_rirs_est, maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5)
)
attr(b_est, "comp_time") <- t_est
attr(b_est, "FPCs") <- mfpca_est
attr(b_est, "nfpc") <- nfpc
save(b_est, file = paste0("simulation/scen_I_110_bamlss_long.Rdata"))



# True FPCs ---------------------------------------------------------------


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

# Load the data
load("simulation/scen_I_130922/data/d110.Rdata")


# Prepare objects for the model on different data sets
mfpca_tru <- MFPCA_cov(cov = cov, basis_funs = b_funs)
mfpca_tru_list <- lapply(seq_along(mfpca_tru$values), 
                         function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
d_rirs_tru <- attach_wfpc(mfpca_tru, d_rirs$data, n = length(mfpca_tru$values))
f_long_tru <- list(
  as.formula(paste0(
    "y ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_along(mfpca_tru$values), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_tru_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker
)


# Model fit
t_est <- system.time(
  b_est <- bamlss(f_long_tru, data = d_rirs_tru, maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5)
)
attr(b_est, "comp_time") <- t_est
attr(b_est, "FPCs") <- mfpca_est
attr(b_est, "nfpc") <- 12
save(b_est, file = paste0("simulation/scen_I_110_bamlss_long_tru1.Rdata"))

