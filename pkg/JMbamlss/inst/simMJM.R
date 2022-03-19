# Data Generation
source("R/simMultiJM.R")
source("R/preprocessing.R")
source("R/eval_mfun.R")
library(tidyverse)

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

dat <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), probmiss = 0.75,
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

aha <- unique(dat$data_short[which(dat$data_short$survtime < 25.1 &
                                     dat$data_short$survtime > 20), "id"])

ggplot(dat$data %>% filter(id %in% aha), aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) + 
  geom_segment(dat$data %>% group_by(id, marker) %>%
                 filter(obstime == max(obstime),
                        id %in% aha),
               mapping = aes(x = obstime, y = y, xend = survtime, yend = y),
               linetype = "dotted") +
  geom_point(dat$data %>% group_by(id, marker) %>%
               filter(obstime == max(obstime),
                      id %in% aha),
             mapping = aes(x = survtime, shape = factor(event)))


mfpca <- preproc_MFPCA(data = dat$data, 
                       uni_mean = "y ~ obstime + x3 + obstime:x3", M = 2, 
                       npc = 2)

data_prep <- attach_wfpc(mfpca, dat$data)


f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, wfpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
f1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, bs = "re") + s(id, obstime, bs = "re"),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
f2 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, marker, bs = "re") + s(id, marker, obstime, bs = "re"),
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
compile_alex()


sink("sim_f1.txt")
b_sim1 <- bamlss(f1, family = mjm_bamlss, data = data_prep, timevar = "obstime",
                verbose_sampler = TRUE)
sink()
save(b_sim1, file = "inst/objects/param_sim01.Rdata")


sink("sim_f.txt")
b_sim <- bamlss(f, family = mjm_bamlss, data = data_prep, timevar = "obstime",
                verbose_sampler = TRUE)
sink()
save(b_sim, file = "inst/objects/param_sim00.Rdata")

sink("sim_f2.txt")
b_sim2 <- bamlss(f2, family = mjm_bamlss, data = data_prep, timevar = "obstime",
                verbose_sampler = TRUE)
sink()
save(b_sim2, file = "inst/objects/param_sim02.Rdata")


# Parametric-like PCA Model -----------------------------------------------

dap <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), probmiss = 0.75,
                  maxfac = 1.5, nmark = 2, param_assoc = FALSE, M = 2, 
                  FPC_bases = NULL, FPC_evals = NULL,
                  mfpc_args = list(type = "split", eFunType = "Poly",
                                   ignoreDeg = NULL, eValType = "linear",
                                   eValScale = 800),
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
                  mu = list(function(t, x){
                    2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                  }, function(t, x, r){
                    2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                  }),
                  sigma = function(t, x) {
                    log(0.6) + 0*t
                  }, 
                  tmax = NULL, seed = 1808, 
                  full = TRUE, file = NULL)

plot(dat$fpc_base)
aha <- unique(dap$data_short[which(dap$data_short$survtime < 25.1 &
                                     dap$data_short$survtime > 20), "id"])

ggplot(dap$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) + 
  geom_segment(dap$data %>% group_by(id, marker) %>%
                 filter(obstime == max(obstime),
                        id %in% aha),
               mapping = aes(x = obstime, y = y, xend = survtime, yend = y),
               linetype = "dotted") +
  geom_point(dap$data %>% group_by(id, marker) %>%
               filter(obstime == max(obstime),
                      id %in% aha),
             mapping = aes(x = survtime, shape = factor(event)))


seq <- seq(0, 25, by = 0.25)
fund <- eFun(argvals = seq, M = 2, type = "Poly")
fpc_base <- multiFunData(
  funData(argvals = seq,
          X = matrix(c(fund@X[1, ],
                       rep(0, length(seq)),
                       fund@X[2, ],
                       rep(0, length(seq))),
                     nrow = 4, byrow = TRUE)),
  funData(argvals = seq,
          X = matrix(c(rep(0, length(seq)),
                       fund@X[1, ],
                       rep(0, length(seq)),
                       fund@X[2, ]),
                     nrow = 4, byrow = TRUE))
)
plot(fpc_base)
norm(fpc_base)

dap <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), probmiss = 0.75,
                  maxfac = 1.5, nmark = 2, param_assoc = FALSE, M = 4, 
                  FPC_bases = fpc_base, FPC_evals = NULL,
                  mfpc_args = list(eValType = "linear", eValScale = 800),
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
                  mu = list(function(t, x){
                    2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                  }, function(t, x, r){
                    2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                  }),
                  sigma = function(t, x) {
                    log(0.6) + 0*t
                  }, 
                  tmax = NULL, seed = 1808, 
                  full = TRUE, file = NULL)

ggplot(dap$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
  