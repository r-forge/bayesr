
# Set Up R Session --------------------------------------------------------

library(MFPCA)
library(tidyverse)
library(Matrix)
library(mvtnorm)
library(survival)
library(bamlss)


# Data Generation
source("R/simMultiJM.R")
source("R/preprocessing.R")

# Helperfunction PCRE
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")
source("R/survint.R")
source("R/compile.R")
source("R/MJM_predict.R")

# Compile the C function
compile_alex()

# Generate data with independent random intercepts
d_indeprirs <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                          probmiss = 0.75, maxfac = 1.5,
                          nmark = 2, param_assoc = TRUE, M = NULL, 
                          FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                          re_cov_mat = matrix(c(0.68, 0, 0, 0,
                                                0, 0.28, 0, 0,
                                                0, 0, 0.68, 0,
                                                0, 0, 0, 0.28), ncol = 4), 
                          ncovar = 2,
                          lambda = function(t, x) {
                            1.37 * t^(0.37)
                          },
                          gamma = function(x) {
                            - 5.8 + 0.48*x[, 3]
                          },
                          alpha = list(function(t, x) {
                            0.64 + 0*t
                          }, function(t, x) {
                            -0.64 + 0*t
                          }),
                          mu = list(function(t, x, r){
                            2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] + 
                              r[, 1] + r[, 2]*t
                          }, function(t, x, r){
                            2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] + 
                              r[, 3] + r[, 4]*t
                          }),
                          sigma = function(t, x) {
                            log(0.6) + 0*t
                          }, 
                          tmax = NULL, seed = 1808, 
                          full = TRUE, file = NULL)


d_indeprirs$data$survtime1 <- d_indeprirs$data$survtime / 25
d_indeprirs$data$obstime1 <- d_indeprirs$data$obstime / 25

set.seed(1808)
seq <- seq(0, 1, by = 0.01)
n <- 100000
b <- mvtnorm::rmvnorm(n = n, mean = c(0,0,0,0), 
                      sigma = diag(c(0.68, 0.28, 0.68, 0.28)))

mfun <- multiFunData(
  funData(argvals = seq,
          X = (b[, 1:2] %*% matrix(c(rep(1, length(seq)), seq),
                                   byrow = TRUE, ncol = length(seq)))),
  funData(argvals = seq,
          X = (b[, 3:4] %*% matrix(c(rep(1, length(seq)), seq),
                                   byrow = TRUE, ncol = length(seq))))
)
mfpca4 <- MFPCA(mFData = mfun, M = 4, 
                uniExpansions = list(list(type = "uFPCA", npc = 2),
                                     list(type = "uFPCA", npc = 2)))

m1 <- list(
  functions = extractObs(mfpca4$functions, 1),
  values = c(1)
)
m2 <- list(
  functions = extractObs(mfpca4$functions, 2),
  values = c(1)
)

dat <- attach_wfpc(mfpca4, d_indeprirs$data, obstime = "obstime1",
                   eval_weight = FALSE)

# Formula
f1 <- list(
  Surv2(survtime1, event, obs = y) ~ -1 + s(survtime1, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime1:marker + x3:marker + obstime1:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre", xt = list("mfpc" = m1)) +
    s(id, fpc.2, bs = "unc_pcre", xt = list("mfpc" = m2)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("indeprirs_unw1.txt")
b_indeprirs_pcre <- bamlss(f1, family = mjm_bamlss, 
                           data = dat, timevar = "obstime1",
                           maxit = 1000, verbose_sampler = TRUE)
sink()
