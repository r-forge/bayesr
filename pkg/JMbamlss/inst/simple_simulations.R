
# Set Up R Session --------------------------------------------------------


library(tidyverse)
library(Matrix)
library(mvtnorm)
library(survival)
library(bamlss)
library(MFPCA)


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

# Compile the C function
compile_alex(dir = "~/Dokumente/Arbeit/Greven/JM/JMbamlss/src/")


# Compare Simple Data Generation Models -----------------------------------

# Generate data with independent random intercepts
d_indepri <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                        probmiss = 0.75, maxfac = 1.5,
                        nmark = 2, param_assoc = TRUE, M = NULL, 
                        FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                        re_cov_mat = matrix(c(0.68, 0, 0, 0.68), ncol = 2), 
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
                            r[, 1]
                        }, function(t, x, r){
                          2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] + 
                            r[, 2]
                        }),
                        sigma = function(t, x) {
                          log(0.6) + 0*t
                        }, 
                        tmax = NULL, seed = 1808, 
                        full = FALSE, file = NULL)

p_indepri <- ggplot(d_indepri, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")


# Appropriate would be to use two independent random effects
f_indepri_onere <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, marker, bs = "re"),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
f_indepri_onere <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, I(marker=="m1"), bs = "re") + s(id, I(marker=="m2"), bs = "re"),
sigma ~ -1 + marker,
alpha ~ -1 + marker
)

# Use an appropriate FPC basis
seq <- seq(0, 25, by = 0.25)
fund <- eFun(argvals = seq, M = 2, type = "Poly")
fpc_base_one <- multiFunData(
  funData(argvals = seq,
          X = matrix(c(fund@X[1, ],
                       rep(0, length(seq))),
                     nrow = 2, byrow = TRUE)),
  funData(argvals = seq,
          X = matrix(c(rep(0, length(seq)),
                       fund@X[1, ]),
                     nrow = 2, byrow = TRUE))
)

mfpca_indepri <- list(
  functions = fpc_base_one,
  values = # Wie sollen die Eigenvalues gewählt werden?
)

# Use different model specifications
f_indepri_onepcre <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, wfpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca_indepri)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

f_indepri_twopcre <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1 bs = "unc_pcre", xt = list("mfpc" = mfpca_indepri)) +
    s(id, wfpc.2 bs = "unc_pcre", xt = list("mfpc" = mfpca_indepri)) ,
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
