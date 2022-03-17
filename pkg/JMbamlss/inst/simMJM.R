# Data Generation
source("R/simMultiJM.R")
source("R/preprocessing.R")
source("R/eval_mfun.R")

dat <- simMultiJM(param_assoc = TRUE, 
                  re_cov_mat = matrix(c(1, 0.5, 0.1, 0.1,
                                        0.5, 1, 0.1, 0.1,
                                        0.1, 0.1, 1, 0.5,
                                        0.1, 0.1, 0.5, 1), ncol = 4),
                  nmark = 2,
                  mu = list(function(time, x, r) {
                    1.25 + r[, 1] + 0.6*sin(x[, 2]) +
                      (-0.01)*time + r[, 2]*time
                  }, function(time, x, r) {
                    -1.25 + r[, 3] + 0.6*sin(x[, 2]) +
                       0.01*time + r[, 4]*time
                  }),
                  full = TRUE)


# Create a scenario as Scenario I in
# Mauff et al. (2020): Joint model with multiple longitudinal outcomes and a
# time-to-event outcome: a corrected two-stage approach

dat <- simMultiJM(nsub = 500, times = seq(0, 25, by = 0.25), probmiss = 0.75,
                  maxfac = 1.5, nmark = 2, param_assoc = TRUE, M = NULL, 
                  FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                  re_cov_mat = matrix(c(0.68, -0.08, 0.1, 0.1, 
                                        -0.08, 0.28, 0.1, 0.1,
                                        0.1, 0.1, 0.68, -0.08,
                                        0.1, 0.1, -0.08, 0.28), ncol = 4), 
                  ncovar = 2,
                  lambda = function(t, x) {
                    1.65 * t^(0.65)
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

mfpca <- preproc_MFPCA(data = dat$data, 
                       uni_mean = "y ~ obstime + x3 + obstime:x3", M = 2, 
                       npc = 2)

data_prep <- attach_wfpc(mfpca, dat$data)


f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, wfpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Helperfunction PCRE
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
compile_alex(dir = "~/Dokumente/Arbeit/Greven/JM/JMbamlss/src/")

b_sim <- bamlss(f, family = mjm_bamlss, data = data_prep, timevar = "obstime")
