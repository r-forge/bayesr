
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

# Compile the C function
compile_alex()


# Compare Simple Data Generation Models -----------------------------------



# Only RI -----------------------------------------------------------------

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
                        full = TRUE, file = NULL)

p_indepri <- ggplot(d_indepri$data, aes(x = obstime, y = y, color = id)) +
  geom_point() +
  geom_line(aes(y = mu)) +
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
f_indepri_byre <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, by = marker, bs = "re"),
sigma ~ -1 + marker,
alpha ~ -1 + marker
)

# Models using random effects
set.seed(1808)
sink("indepri_onere.txt")
b_indepri_onere <- bamlss(f_indepri_onere, family = mjm_bamlss, 
                          data = d_indepri$data, timevar = "obstime",
                          verbose_sampler = TRUE)
sink()

set.seed(1808)
sink("indepri_byre.txt")
b_indepri_byre <- bamlss(f_indepri_byre, family = mjm_bamlss, 
                          data = d_indepri$data, timevar = "obstime",
                          verbose_sampler = TRUE)
sink()


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
  values = c(17, 17)
)


# Attach the FPCS
d_indepri$data <- attach_wfpc(mfpca_indepri, d_indepri$data)

# Use principal components
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
    s(id, wfpc.1, bs = "unc_pcre",
      xt = list("mfpc" = list(functions = extractObs(fpc_base_one, 1),
                              values = c(17)))) +
    s(id, wfpc.2, bs = "unc_pcre",
      xt = list("mfpc" = list(functions = extractObs(fpc_base_one, 2),
                              values = c(17)))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Models using PCREs
set.seed(1808)
sink("indepri_onepcre.txt")
b_indepri_onepcre <- bamlss(f_indepri_onepcre, family = mjm_bamlss, 
                            data = d_indepri$data, timevar = "obstime",
                            verbose_sampler = TRUE)
sink()

set.seed(1808)
sink("indepri_twopcre.txt")
b_indepri_twopcre <- bamlss(f_indepri_twopcre, family = mjm_bamlss, 
                            data = d_indepri$data, timevar = "obstime",
                            verbose_sampler = TRUE)
sink()

save(b_indepri_onere, b_indepri_byre, b_indepri_onepcre, b_indepri_twopcre,
     file = "inst/objects/indepri_models.Rdata")

# Compare models
load("inst/objects/indepri_models.Rdata")
set.seed(1808)
ids <- sample(unique(d_indepri$data_short[which(
  d_indepri$data_short$survtime < 25.1 &d_indepri$data_short$survtime > 20), 
  "id"]), 5, replace = FALSE)

p_re <- ggplot(data = d_indepri$data %>% 
                 mutate(fit_onere = b_indepri_onere$fitted.values$mu, 
                        fit_byre = b_indepri_byre$fitted.values$mu) %>% 
                 filter(id %in% ids), aes(x = obstime, color = id)) +
  geom_point(aes(y = y), size = 1.3) +
  geom_line(aes(y = fit_onere), linetype = "dotted") +
  geom_line(aes(y = fit_byre), linetype = "dashed") +
  geom_line(aes(y = mu)) +
  facet_grid(~marker)
p_pcre <- ggplot(data = d_indepri$data %>% 
                   mutate(fit_onepcre = b_indepri_onepcre$fitted.values$mu, 
                          fit_twopcre = b_indepri_twopcre$fitted.values$mu) %>% 
                   filter(id %in% ids), aes(x = obstime, color = id)) +
  geom_point(aes(y = y), size = 1.3) +
  geom_line(aes(y = fit_onepcre), linetype = "dotted") +
  geom_line(aes(y = fit_twopcre), linetype = "dashed") +
  geom_line(aes(y = mu)) +
  facet_grid(~marker)

# Models are actually equal
all.equal(b_indepri_byre$fitted.values$mu, b_indepri_onere$fitted.values$mu)
all.equal(b_indepri_onere$fitted.values$mu, b_indepri_onepcre$fitted.values$mu)
all.equal(b_indepri_byre$fitted.values$mu, b_indepri_twopcre$fitted.values$mu)

summary(b_indepri_byre$samples[[
  1]][, grep("accepted", colnames(b_indepri_byre$samples[[1]]))])$statistics
summary(b_indepri_onere$samples[[
  1]][, grep("accepted", colnames(b_indepri_onere$samples[[1]]))])$statistics
summary(b_indepri_twopcre$samples[[
  1]][, grep("accepted", colnames(b_indepri_twopcre$samples[[1]]))])$statistics
summary(b_indepri_onepcre$samples[[
  1]][, grep("accepted", colnames(b_indepri_onepcre$samples[[1]]))])$statistics



# Independent RI + RS -----------------------------------------------------

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

p_indeprirs <- ggplot(d_indeprirs$data, aes(x = obstime, y = y, color = id)) +
  geom_point() +
  geom_line(aes(y = mu)) +
  facet_grid(~marker) +
  theme(legend.position = "none")

