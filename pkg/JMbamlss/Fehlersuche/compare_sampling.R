
# Multivariate Joint Model Example ----------------------------------------

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

# Alternative bamlss Code (more options)
source("Fehlersuche/JM.R")
library(Matrix)
library(mvtnorm)
library(survival)
library(bamlss)
library(MFPCA)

load("inst/objects/simu_meike/a_linear_constant_zero_150_1.RData")

f_jm <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, k = 5, bs = "ps"),
  gamma ~ s(x1, k = 5, bs = "ps"),
  mu ~ obstime + s(x2, k = 5, bs = "ps") + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_jm <- bamlss(f_jm, data = d$data, family = "jmFEHLERSUCHE", 
               timevar = "obstime", idvar = "id", subdivisions = 25,
               maxit = 1500)

b_jm_lambda <- bamlss(f_jm, data = d$data, family = "jmFEHLERSUCHE", 
                      timevar = "obstime", idvar = "id", subdivisions = 25,
                      optimizer = FALSE, start = parameters(b_jm), 
                      prop_pred = "lambda")
b_jm_gamma <- bamlss(f_jm, data = d$data, family = "jmFEHLERSUCHE", 
                     timevar = "obstime", idvar = "id", subdivisions = 25,
                     optimizer = FALSE, start = parameters(b_jm), 
                     prop_pred = "gamma")

f_mjm <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, k = 5, bs = "ps"),
  gamma ~ s(x1, k = 5, bs = "ps"),
  mu ~ obstime + s(x2, k = 5, bs = "ps") + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)

b_mjm <- bamlss(f_mjm, data = d$data, family = "mjm_bamlss", 
                timevar = "obstime", optimizer = FALSE,
                start = parameters(b_jm))

b_mjm_lambda <- bamlss(f_mjm, data = d$data, family = "mjm_bamlss", 
                       timevar = "obstime", optimizer = FALSE,
                       start = parameters(b_jm), 
                       prop_pred = "lambda")

b_mjm_gamma <- bamlss(f_mjm, data = d$data, family = "mjm_bamlss", 
                      timevar = "obstime", optimizer = FALSE,
                      start = parameters(b_jm), 
                      prop_pred = "gamma")

b_mjm_lgs <- bamlss(f_mjm, data = d$data, family = "mjm_bamlss", 
                    timevar = "obstime", optimizer = FALSE,
                    start = parameters(b_jm), 
                    prop_pred = c("lambda", "gamma", "sigma"))

# Compare accepted ratios
summary(b_jm$samples[[1]][, c(8, 16, 20, 28, 182, 187, 191, 195)])
summary(b_mjm$samples[[1]][, c(8, 16, 20, 28, 182, 187, 191, 195)])
summary(b_jm_lambda$samples[[1]][, 8])
summary(b_mjm_lambda$samples[[1]][, 8])
summary(b_jm_gamma$samples[[1]][, c(16, 20)])
summary(b_mjm_gamma$samples[[1]][, c(8, 12)])
summary(b_mjm_lgs$samples[[1]][, c(8, 16, 20, 24)])
