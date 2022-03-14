
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

# Data Generation ---------------------------------------------------------

if(!exists("d")) {
  
  library("survival")
  library("bamlss")
  library("MFPCA")
  
  
  set.seed(1808)
  d <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121),
                  probmiss = 0.75, maxfac = 1.5, nmark = 2, M = 6, ncovar = 2,
                  lambda = function(t, x) {
                    1.4*log((120*t + 10)/1000)
                  },
                  alpha = rep(list(function(t, x) {
                    0.3 + 0*t
                  }), 2),
                  mu = rep(list(function(t, x){
                    1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                  }), 2),
                  mfpc_args = list(type = "split", eFunType = "Poly",
                                   ignoreDeg = NULL, eValType = "linear",
                                   eValScale = 1))
  mfpca <- preproc_MFPCA(data = d, uni_mean = "y ~ s(obstime) + s(x2)", M = 2)
  d_simp <- simMultiJM(nsub = 300, times = seq(0, 1, length.out = 121),
                       probmiss = 0.75, maxfac = 1.5, nmark = 1, M = 6, 
                       ncovar = 2,
                       lambda = function(t, x) {
                         1.4*log((120*t + 10)/1000)
                       },
                       alpha = list(function(t, x) {
                         0.3 + 0*t
                       }),
                       mu = list(function(t, x){
                         1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                       }),
                       sigma = function(t, x) {-50 + 0*t}, 
                       full = TRUE,
                       mfpc_args = list(type = "split", eFunType = "Poly",
                                        ignoreDeg = NULL, eValType = "linear",
                                        eValScale = 1))
  
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

sink("MJM_pcre.txt")
b_sample <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime")
sink()

## Problem bei Alpha-Prädiktor:
# Warum sollte sich die Likelihood so extrem ändern, nur weil andere alpha-
# Parameter verwendet werden? Zumal die Koeffizienten jetzt nicht heillos unter-
# schiedlich sind. 
# Wird vielleicht doch einer der longitudinalen Prädiktoren falsch upgedated?
# Allerdings, warum sollte es dann zwischendurch immer mal wieder funktionieren?



# RE Model ----------------------------------------------------------------

f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + ti(obstime, by = marker) +
    ti(id, bs = "re") +
    ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
sink("MJM_re.txt")
b <- bamlss(f_re, family = mjm_bamlss, data = d, timevar = "obstime")
sink()


# Very Simple Example ----------------------- ------------------------------

# set.seed(1808)
# d_simp <- simMultiJM(nsub = 250, times = seq(0, 1, length.out = 121),
#                      nmark = 1, full = TRUE,
#                      lambda = function(t, x) {
#                         1.4*log((120*t + 10)/1000)
#                      },
#                      alpha = list(function(t, x) {
#                         0.3 + 0*t
#                      }),
#                      mu = list(function(t, x){
#                         1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
#                      }),
#                      sigma = function(t, x) {-50 + 0*t})
# ggplot(d_simp$data, aes(x = obstime, y = y, group = id)) + 
#   geom_line() +
#   geom_segment(d_simp$data %>% group_by(id) %>% 
#                  filter(obstime == max(obstime)),
#                mapping = aes(x = obstime, y = y, xend = survtime, yend = y),
#                linetype = "dotted") +
#   geom_point(d_simp$data %>% group_by(id) %>% 
#                filter(obstime == max(obstime)),
#              mapping = aes(x = survtime, shape = factor(event)))

f_simp <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ s(obstime) + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)

# b_simp <- bamlss(f_simp, family = mjm_bamlss, data = d_simp$data,
#                  timevar = "obstime", sampler = FALSE)
# load("inst/objects/m_simp_opt.Rdata")
b_simp_samp <- bamlss(f_simp, family = mjm_bamlss, data = d_simp$data, 
                 timevar = "obstime", optimizer = FALSE, 
                 start = parameters(b_simp_jm))

f_simp_jm <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ s(obstime) + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
# 
debug(bamlss::sam_JM)
b_simp_jm <- bamlss(f_simp_jm , family = "jm", data = d_simp$data,
                    timevar = "obstime", idvar = "id", prop_pred = "lambda")
