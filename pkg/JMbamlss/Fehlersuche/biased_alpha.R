
# Biased Alpha Estimation -------------------------------------------------


# Data Generation and Preprocessing
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

# Alternative bamlss Code (more options)
#source("Fehlersuche/JM.R")

library(Matrix)
library(tidyverse)
library(refund)
library(survival)
library(bamlss)
library(funData)


# Linear BH, Small Sigma --------------------------------------------------

set.seed(1808)

# Linear hazard, no association, constant marker values (uncentered)
dat_unc <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                         lambda = function(t, x) -5 + 0.02*t,
                         gamma = function(x) 0,
                         alpha = list(function(t, x) 0*t),
                         mu = list(function(t, x) 1.25),
                         sigma = function(t, x) -50 + 0*t,
                         full = TRUE)
ggplot(dat_unc$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))



# Bamlss ------------------------------------------------------------------

f_bamlss <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

m_bamlss <- bamlss(f_bamlss, family = "jm", data = dat_unc$data, 
                   timevar = "obstime", idvar = "id", sampler = FALSE)
m_bamlss$parameters$alpha



# Own Implementation ------------------------------------------------------

f_mjm <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
m_mjm <- bamlss(f_mjm, family = "mjm", data = dat_unc$data, 
                   timevar = "obstime", sampler = FALSE)
m_mjm$parameters$alpha

