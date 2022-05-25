
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


# Independent RI + RS -----------------------------------------------------

# Generate data with independent random intercepts
d_25 <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
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
d_1 <- d_25$data %>% mutate(survtime = survtime / 25, obstime = obstime / 25)


# Generate true MFPCs
seq25 <- seq(0, 25, by = 0.25)
seq1 <- seq(0, 1, by = 0.01)
set.seed(1808)
n <- 100000
b <- mvtnorm::rmvnorm(n = n, mean = c(0,0,0,0), 
                      sigma = diag(c(0.68, 0.28, 0.68, 0.28)))

mfun25 <- multiFunData(
  funData(argvals = seq25,
          X = (b[, 1:2] %*% matrix(c(rep(1, length(seq25)), seq25),
                                   byrow = TRUE, ncol = length(seq25)))),
  funData(argvals = seq25,
          X = (b[, 3:4] %*% matrix(c(rep(1, length(seq25)), seq25),
                                   byrow = TRUE, ncol = length(seq25))))
)

mfun1 <- multiFunData(
  funData(argvals = seq1,
          X = (b[, 1:2] %*% matrix(c(rep(1, length(seq25)), seq25),
                                   byrow = TRUE, ncol = length(seq25)))),
  funData(argvals = seq1,
          X = (b[, 3:4] %*% matrix(c(rep(1, length(seq25)), seq25),
                                   byrow = TRUE, ncol = length(seq25))))
)


mfpca25_2 <- MFPCA(mFData = mfun25, M = 4, 
                       uniExpansions = list(list(type = "uFPCA"),
                                            list(type = "uFPCA")))
mfpca25_4 <- MFPCA(mFData = mfun25, M = 4, 
                       uniExpansions = list(list(type = "uFPCA", npc = 2),
                                            list(type = "uFPCA", npc = 2)))
mfpca1_2 <- MFPCA(mFData = mfun1, M = 4,
                      uniExpansions = list(list(type = "uFPCA"),
                                           list(type = "uFPCA")))
mfpca1_4 <- MFPCA(mFData = mfun1, M = 4,
                      uniExpansions = list(list(type = "uFPCA", npc = 2),
                                           list(type = "uFPCA", npc = 2)))

mfpca25_2_list <- lapply(1:2, function (i, mfpca = mfpca25_2) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
mfpca25_4_list <- lapply(1:4, function (i, mfpca = mfpca25_4) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
mfpca1_2_list <- lapply(1:2, function (i, mfpca = mfpca1_2) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
mfpca1_4_list <- lapply(1:4, function (i, mfpca = mfpca1_4) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})



d_25_w4 <- attach_wfpc(mfpca25_4, d_25$data, eval_weight = TRUE)
d_25_4 <- attach_wfpc(mfpca25_4, d_25$data, eval_weight = FALSE)
d_1_w4 <- attach_wfpc(mfpca1_4, d_1, eval_weight = TRUE)
d_1_4 <- attach_wfpc(mfpca1_4, d_1, eval_weight = FALSE)

d_25_w2 <- attach_wfpc(mfpca25_2, d_25$data, eval_weight = TRUE)
d_25_2 <- attach_wfpc(mfpca25_2, d_25$data, eval_weight = FALSE)
d_1_w2 <- attach_wfpc(mfpca1_2, d_1, eval_weight = TRUE)
d_1_2 <- attach_wfpc(mfpca1_2, d_1, eval_weight = FALSE)

# Different Starting Model Specifications ---------------------------------


# Starting model
f_wpc4_s1_25 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, wfpc.2, wfpc.3, wfpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4, "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc4_s1_25.txt")
b_wpc4_s1_25 <- bamlss(f_wpc4_s1_25, family = mjm_bamlss, 
                       data = d_25_w4, 
                       timevar = "obstime",
                       maxit = 1000, verbose_sampler = TRUE)
sink()

# Starting model without weighting
f_pc4_s1_25 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, fpc.2, fpc.3, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4, "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc4_s1_25.txt")
b_pc4_s1_25 <- bamlss(f_pc4_s1_25, family = mjm_bamlss, 
                      data = d_25_4, 
                      timevar = "obstime",
                      maxit = 1000, verbose_sampler = TRUE)
sink()

# Starting model with 4 variance parameters
f_wpc4_s4_25 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[1]], "eval_weight" = TRUE)) +
    s(id, wfpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[2]], "eval_weight" = TRUE)) +
    s(id, wfpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[3]], "eval_weight" = TRUE)) +
    s(id, wfpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[4]], "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc4_s4_25.txt")
b_wpc4_s4_25 <- bamlss(f_wpc4_s4_25, family = mjm_bamlss, 
                       data = d_25_w4, 
                       timevar = "obstime",
                       maxit = 1000, verbose_sampler = TRUE)
sink()

# Starting model with interval 1
f_wpc4_s1_1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, wfpc.2, wfpc.3, wfpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4, "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc4_s1_1.txt")
b_wpc4_s1_1 <- bamlss(f_wpc4_s1_1, family = mjm_bamlss, 
                       data = d_1_w4, 
                       timevar = "obstime",
                       maxit = 1000, verbose_sampler = TRUE)
sink()



# Different Models with Separate S ----------------------------------------

# Unweighted FPCs
f_pc4_s4_25 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[2]], "eval_weight" = FALSE)) +
    s(id, fpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[3]], "eval_weight" = FALSE)) +
    s(id, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_list[[4]], "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc4_s4_25.txt")
b_pc4_s4_25 <- bamlss(f_pc4_s4_25, family = mjm_bamlss, 
                       data = d_25_4, 
                       timevar = "obstime",
                       maxit = 1000, verbose_sampler = TRUE)
sink()

# Using only 2 FPCs
f_wpc2_s2_25 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_2_list[[1]], "eval_weight" = TRUE)) +
    s(id, wfpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_2_list[[2]], "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc2_s2_25.txt")
b_wpc2_s2_25 <- bamlss(f_wpc2_s2_25, family = mjm_bamlss, 
                       data = d_25_w2, 
                       timevar = "obstime",
                       maxit = 1000, verbose_sampler = TRUE)
sink()

# On Scale to 1
f_wpc4_s4_1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, wfpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[1]], "eval_weight" = TRUE)) +
    s(id, wfpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[2]], "eval_weight" = TRUE)) +
    s(id, wfpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[3]], "eval_weight" = TRUE)) +
    s(id, wfpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[4]], "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc4_s4_1.txt")
b_wpc4_s4_1 <- bamlss(f_wpc4_s4_1, family = mjm_bamlss, 
                      data = d_1_w4, 
                      timevar = "obstime",
                      maxit = 1000, verbose_sampler = TRUE)
sink()


# Different Models with All Combinations ----------------------------------

# On Scale to 1
f_pc4_s4_1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[2]], "eval_weight" = FALSE)) +
    s(id, fpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[3]], "eval_weight" = FALSE)) +
    s(id, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[4]], "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc4_s4_1.txt")
b_pc4_s4_1 <- bamlss(f_pc4_s4_1, family = mjm_bamlss, 
                      data = d_1_4, 
                      timevar = "obstime",
                      maxit = 1000, verbose_sampler = TRUE)
sink()


# Only two FPCs
f_pc2_s2_1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_2_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_2_list[[2]], "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc2_s2_1.txt")
b_pc2_s2_1 <- bamlss(f_pc2_s2_1, family = mjm_bamlss, 
                     data = d_1_2, 
                     timevar = "obstime",
                     maxit = 1000, verbose_sampler = TRUE)
sink()



# Different Models --------------------------------------------------------


# Using only 2 FPCs
f_wpc2_s2_25_RI <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, bs = "re", by = marker) +
    s(id, wfpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_2_list[[1]], "eval_weight" = TRUE)) +
    s(id, wfpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_2_list[[2]], "eval_weight" = TRUE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_wpc2_s2_25_RI.txt")
b_wpc2_s2_25_RI <- bamlss(f_wpc2_s2_25_RI, family = mjm_bamlss, 
                          data = d_25_w2, 
                          timevar = "obstime",
                          maxit = 1000, verbose_sampler = TRUE)
sink()

# Model using manual FPC basis
fund25 <- eFun(argvals = seq25, M = 2, type = "Poly")
fpc_base_two <- multiFunData(
  funData(argvals = seq25,
          X = matrix(c(fund25@X[1, ],
                       rep(0, length(seq25)),
                       fund25@X[2, ],
                       rep(0, length(seq25))),
                     nrow = 4, byrow = TRUE)),
  funData(argvals = seq25,
          X = matrix(c(rep(0, length(seq25)),
                       fund25@X[1, ],
                       rep(0, length(seq25)),
                       fund25@X[2, ]),
                     nrow = 4, byrow = TRUE))
)
mfpca25_4_man <- list(
  functions = fpc_base_two,
  values = c(1, 1, 1, 1)
)
mfpca25_4_man_list <- lapply(1:4, function (i, mfpca = mfpca25_4_man) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_pc4_s4_25_MAN <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_man_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_man_list[[2]], "eval_weight" = FALSE)) +
    s(id, fpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_man_list[[3]], "eval_weight" = FALSE)) +
    s(id, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca25_4_man_list[[4]], "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
d_25_4_MAN <- attach_wfpc(mfpca25_4_man, d_25$data, eval_weight = FALSE)

set.seed(1808)
sink("rirs/rirs_pc4_s4_25_MAN.txt")
b_pc4_s4_25_MAN <- bamlss(f_pc4_s4_25_MAN, family = mjm_bamlss, 
                      data = d_25_4_MAN, 
                      timevar = "obstime",
                      maxit = 1000, verbose_sampler = TRUE)
sink()


# Model using additional RI
f_pc4_s4_1_RE <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[2]], "eval_weight" = FALSE)) +
    s(id, fpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[3]], "eval_weight" = FALSE)) +
    s(id, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca1_4_list[[4]], "eval_weight" = FALSE)) +
    s(id, bs = "re", by = marker),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc4_s4_1_RE.txt")
b_pc4_s4_1_RE <- bamlss(f_pc4_s4_1_RE, family = mjm_bamlss, 
                        data = d_1_4, 
                        timevar = "obstime",
                        maxit = 1000, verbose_sampler = TRUE)
sink()

# Model using only RI+RS
f_rirs <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, by = marker, bs = "re") + s(id, obstime, by =  marker, bs = "re"),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_re.txt")
b_rirs <- bamlss(f_rirs, family = mjm_bamlss, 
                 data = d_1_4, timevar = "obstime",
                 maxit = 1000, verbose_sampler = TRUE)
sink()


# Model with estimated FPCs -----------------------------------------------

mfpca_es <- preproc_MFPCA(d_1, uni_mean = "y ~ 1+obstime + x3 + obstime:x3",
                          M = 4, npc = 2)
mfpca_es_list <- lapply(1:4, function (i, mfpca = mfpca_es) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
d_1_es <- attach_wfpc(mfpca_es, d_1, eval_weight = FALSE)


# On Scale to 1
f_pc4_s4_1_es <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
    s(id, fpc.1, bs = "unc_pcre",
      xt = list("mfpc" = mfpca_es_list[[1]], "eval_weight" = FALSE)) +
    s(id, fpc.2, bs = "unc_pcre",
      xt = list("mfpc" = mfpca_es_list[[2]], "eval_weight" = FALSE)) +
    s(id, fpc.3, bs = "unc_pcre",
      xt = list("mfpc" = mfpca_es_list[[3]], "eval_weight" = FALSE)) +
    s(id, fpc.4, bs = "unc_pcre",
      xt = list("mfpc" = mfpca_es_list[[4]], "eval_weight" = FALSE)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1808)
sink("rirs/rirs_pc4_s4_1_es.txt")
b_pc4_s4_1_es <- bamlss(f_pc4_s4_1_es, family = mjm_bamlss, 
                        data = d_1_es,
                        timevar = "obstime",
                        maxit = 1000, verbose_sampler = TRUE)
sink()



# Model using different Data Set ------------------------------------------

# Generate data with independent random intercepts
d_25_other <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
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
                         tmax = NULL, seed = 1005, 
                         full = TRUE, file = NULL)
d_1_other <- d_25_other$data %>%
  mutate(survtime = survtime / 25, obstime = obstime / 25)

d_1_4 <- attach_wfpc(mfpca1_4, d_1_other, eval_weight = FALSE)
set.seed(1808)
sink("rirs/rirs_pc4_s4_1_other.txt")
b_pc4_s4_1_other <- bamlss(f_pc4_s4_1, family = mjm_bamlss, 
                           data = d_1_4_other, 
                           timevar = "obstime",
                           maxit = 1000, verbose_sampler = TRUE)
sink()


# MGCV Modell -------------------------------------------------------------

g <- gam(y ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
           s(id, wfpc.1, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[1]], "eval_weight" = TRUE)) +
           s(id, wfpc.2, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[2]], "eval_weight" = TRUE)) +
           s(id, wfpc.3, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[3]], "eval_weight" = TRUE)) +
           s(id, wfpc.4, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[4]], "eval_weight" = TRUE)), 
         data = d_25_w4)

g2 <- gam(y ~ -1 + marker + obstime:marker + x3:marker + obstime:marker:x3 +
           s(id, wfpc.1, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[1]], "eval_weight" = TRUE)) +
           s(id, wfpc.2, bs = "unc_pcre",
             xt = list("mfpc" = mfpca25_4_list[[2]], "eval_weight" = TRUE)), 
         data = d_25_w4)

ggplot(d_25$data %>% mutate(fit4 = g$fitted.values),
       aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = fit4), alpha = 0.7) +
  geom_point(aes(y = y)) +
  facet_wrap(~ marker, scales = "free") +
  theme(legend.position = "none")

ggplot(d_25$data %>% mutate(fit2 = g2$fitted.values),
       aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu)) +
  geom_line(aes(y = fit2), alpha = 0.7) +
  geom_point(aes(y = y)) +
  facet_wrap(~ marker, scales = "free") +
  theme(legend.position = "none")

# Fits are reasonable, but the two-pcre fit is noticeably worse

 

# Does the smooth construct use additional centering? ---------------------

debug(MJM_mcmc)
b_wpc4_s1_25 <- bamlss(f_wpc4_s1_25, family = mjm_bamlss, 
                       data = d_25_w4, 
                       timevar = "obstime", optimizer = FALSE, 
                       start = parameters(b_wpc4_s1_25))

