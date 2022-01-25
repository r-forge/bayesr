
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


# Data Generation ---------------------------------------------------------

if(!exists("d")) {
  
  library("survival")
  library("bamlss")
  library("MFPCA")
  
  
  set.seed(1808)
  d <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121), 
                  lambda = function(t, x) {
                    1.4*log((120*t + 10)/1000)
                  },
                  alpha = rep(list(function(t, x) {
                    0.3 + 0*t
                  }), 2),
                  mu = rep(list(function(t, x){
                    1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                  }), 2))
  mfpca <- preproc_MFPCA(data = d, uni_mean = "y ~ s(obstime) + s(x2)", M = 2)
}




# PCRE Model --------------------------------------------------------------

f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

b <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime",
            sampler = FALSE, maxit = 500)
# load("inst/objects/m_optim.Rdata")

b_sample <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime",
                   optimizer = FALSE, start = parameters(b))
#, n.iter = 4, burnin = 1, step = 1)

## Problem bei Alpha-Prädiktor:
# Warum sollte sich die Likelihood so extrem ändern, nur weil andere alpha-
# Parameter verwendet werden? Zumal die Koeffizienten jetzt nicht heillos unter-
# schiedlich sind. 
# Wird vielleicht doch einer der longitudinalen Prädiktoren falsch upgedated?
# Allerdings, warum sollte es dann zwischendurch immer mal wieder funktionieren?



# RE Model ----------------------------------------------------------------

# f_re <- list(
#   Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
#   gamma ~ 1, 
#   mu ~ -1 + marker + ti(obstime, by = marker) +
#     ti(id, bs = "re") +
#     ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
#   sigma ~ -1 + marker,
#   alpha ~ -1 + marker + s(survtime, by = marker)
# )
# b <- bamlss(f_re, family = mjm_bamlss, data = d, timevar = "obstime")


