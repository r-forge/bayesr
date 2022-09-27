
# Set Up Small Simulation -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


# Always
library(survival)
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/MJM_opt.R")
source("R/opt_updating_cpp.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing_cpp.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/fpca.R")
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")


# Setting for the simulation
start <- 100
stop <- 119
number_cores <- 2
setting <- "scen_I"
dir.create(paste0("simulation/", setting), showWarnings = FALSE)
dir.create(paste0("simulation/", setting, "/data"), showWarnings = FALSE)
dir.create(paste0("simulation/", setting, "/bamlss_tru"), showWarnings = FALSE)
dir.create(paste0("simulation/", setting, "/bamlss_est"), showWarnings = FALSE)
dir.create(paste0("simulation/", setting, "/jmb"), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Objects for all clusters ------------------------------------------------

# Number of individuals
n <- 150

# Covariance matrix for the data generation
auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
cross <- matrix(rep(0.03, 4), ncol = 2)
cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
                1, 0, 1, 0.75, 0.5, 0,
                0.75, 1, 0, 1, 0.75, 0.5,
                0.5, 0.75, 1, 0, 1, 0.75,
                0, 0.5, 0.75, 1, 0, 1,
                0, 0, 0.5, 0.75, 1, 0),
              ncol = 6)
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)), auto)

# Use true dependence structure to get a good estimate of the underlying FPC
# basis
seq1 <- seq(0, 1, by = 0.01)
set.seed(1808)
n_pca <- 100000
b <- mvtnorm::rmvnorm(n = n_pca, sigma = cov)
mfun1 <- multiFunData(
  funData(argvals = seq1,
          X = (b[, 1:2] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1)))),
  funData(argvals = seq1,
          X = (b[, 3:4] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1)))),
  funData(argvals = seq1,
          X = (b[, 5:6] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1)))),
  funData(argvals = seq1,
          X = (b[, 7:8] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1)))),
  funData(argvals = seq1,
          X = (b[, 9:10] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1)))),
  funData(argvals = seq1,
          X = (b[, 11:12] %*% matrix(c(rep(1, length(seq1)), seq1),
                                   byrow = TRUE, ncol = length(seq1))))
  )
mfpca_tru <- MFPCA(mFData = mfun1, M = 12,
                   uniExpansions = list(list(type = "uFPCA", npc = 2),
                                        list(type = "uFPCA", npc = 2),
                                        list(type = "uFPCA", npc = 2),
                                        list(type = "uFPCA", npc = 2),
                                        list(type = "uFPCA", npc = 2),
                                        list(type = "uFPCA", npc = 2)))
save(mfpca_tru, file = paste0("simulation/", setting, "/mfpca_tru.Rdata"))

# Prepare the model using the true FPCs
mfpca <- lapply(1:7, function (i, m = mfpca_tru) {
  list(functions = extractObs(m$functions, i),
       values = m$values[i])
})
f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +
    s(id, fpc.1, bs = "unc_pcre", xt = list("mfpc" = mfpca[[1]])) + 
    s(id, fpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca[[2]])) + 
    s(id, fpc.3, bs = "unc_pcre", xt = list("mfpc" = mfpca[[3]])) +
    s(id, fpc.4, bs = "unc_pcre", xt = list("mfpc" = mfpca[[4]])) + 
    s(id, fpc.5, bs = "unc_pcre", xt = list("mfpc" = mfpca[[5]])) + 
    s(id, fpc.6, bs = "unc_pcre", xt = list("mfpc" = mfpca[[6]])) +
    s(id, fpc.7, bs = "unc_pcre", xt = list("mfpc" = mfpca[[7]])),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)



# Simulation function -----------------------------------------------------

parallel_mods <- function(i) {
  set.seed(i)
  
  # Simulate the data
  d_rirs <- simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                       probmiss = 0.75, maxfac = 1.75,
                       nmark = 6, long_assoc = "param", M = NULL, 
                       FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                       re_cov_mat = cov,
                       ncovar = 2,
                       lambda = function(t, x) {
                         1.37 * t^(0.37)
                       },
                       gamma = function(x) {
                         -1.5 + 0.48*x[, 3]
                       },
                       alpha = list(function(t, x) {
                         1.5 + 0*t
                       }, function(t, x) {
                         0.6 + 0*t
                       }, function(t, x) {
                         0.3 + 0*t
                       }, function(t, x) {
                         -0.3 + 0*t
                       }, function(t, x) {
                         -0.6 + 0*t
                       }, function(t, x) {
                         -1.5 + 0*t
                       }),
                       mu = list(function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 1] + r[, 2]*t
                       }, function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 3] + r[, 4]*t
                       }, function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 5] + r[, 6]*t
                       }, function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 7] + r[, 8]*t
                       }, function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 9] + r[, 10]*t
                       }, function(t, x, r){
                         0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                           r[, 11] + r[, 12]*t
                       }),
                       sigma = function(t, x) {
                         log(0.06) + 0*t
                       }, 
                       tmax = NULL, seed = NULL, 
                       full = TRUE, file = NULL)
  
  save(d_rirs, file = paste0("simulation/", setting, "/data/d", i, ".Rdata"))
  # ggplot(data = d_rirs$data_short %>% slice_head(n = n),
  #        aes(x = id, y = survtime, colour = as.factor(event))) +
  #   geom_point()
  # ggplot(d_rirs$data, aes(x = obstime, y = y, color = id)) +
  #   geom_line() +
  #   facet_grid(~marker) +
  #   theme(legend.position = "none")
  # ggplot(d_rirs$data_hypo %>% filter(event == 1), aes(x = obstime, y = mu, color = id)) +
  #   geom_line() +
  #   facet_grid(~marker) +
  #   theme(legend.position = "none")
  # summary(c(table(d_rirs$data$id)))
  # table(d_rirs$data_short[seq_len(n), "event"])
  
  # Prepare the model using estimated FPCs
  # mfpca_es <- preproc_MFPCA(d_rirs$data, 
  #                           uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
  #                           method = "FPCA")
  
  # Estimate the model using the true FPCs
  # d_rirs_tru <- subset(d_rirs$data, 
  #                      select = -grep("fpc\\.", colnames(d_rirs$data)))
  d_rirs_tru <- attach_wfpc(mfpca_tru, d_rirs$data, n = 7)
  
  t_tru <- system.time(
    b_tru <- try(bamlss(f_tru, family = mjm_bamlss, data = d_rirs_tru, 
                        timevar = "obstime", maxit = 1500, n.iter = 5500,
                        burnin = 500, thin = 5)))
  attr(b_tru, "comp_time") <- t_tru
  save(b_tru, file = paste0("simulation/", setting, "/bamlss_tru/b", i, 
                            ".Rdata"))
  
  # Estimate the model using estimated FPCs
  try({
    few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                     function (x) any(x < 2))
    mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                                filter(id %in% paste(which(!few_obs))) %>% 
                                droplevels(), 
                              uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                              npc = 2, nbasis = 4)
    vals <- which(mfpca_est$values > 0)
    nfpc <- min(which(
      cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
    mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
      list(functions = extractObs(mfpca$functions, i),
           values = mfpca$values[i])
    })
    # d_rirs_est <- subset(d_rirs$data, 
    #                     select = -grep("fpc\\.", colnames(d_rirs$data)))
    d_rirs_est <- attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
    f_est <- list(
      Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
      gamma ~ 1 + x3,
      as.formula(paste0(
        "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
        paste0(lapply(seq_len(nfpc), function(x) {
          paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
                 "mfpca_est_list[[", x, "]]))")
        }), collapse = " + "))),
      sigma ~ -1 + marker,
      alpha ~ -1 + marker
    )
    t_est <- system.time(
      b_est <- bamlss(f_est, family = mjm_bamlss, data = d_rirs_est, 
                      timevar = "obstime", maxit = 1500, n.iter = 5500,
                      burnin = 500, thin = 5)
    )
    attr(b_est, "comp_time") <- t_est
    
  })
  save(b_est, file = paste0("simulation/", setting, "/bamlss_est/b", i, 
                            ".Rdata"))

  
  # Estimate the model using JMbayes
  d_rirs_jmb <- d_rirs$data %>% 
    pivot_wider(names_from = marker, values_from = y)
  CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                  data = d_rirs$data_short %>% filter(marker == "m1"))
  lm1 <- lme(m1 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)
  lm2 <- lme(m2 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)
  lm3 <- lme(m3 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)
  lm4 <- lme(m4 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)
  lm5 <- lme(m5 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)
  lm6 <- lme(m6 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit)

  t_jm <- system.time(
    jmb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "obstime",
              n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
              control = list(cores = 1))
  )
  attr(jmb, "comp_time") <- t_jm
  save(jmb, file = paste0("simulation/", setting, "/jmb/jmb_", i, 
                            ".Rdata"))
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_mods, mc.cores = number_cores)
# 
# # Actual parallel section
# cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
# cl <- makeForkCluster(number_cores)
# sim <- parLapply(cl, start:stop, simulate_models)
# stopCluster(cl)
# Sys.time()
# 
# rm(list =ls())
# quit("no")

# 
# # Try scenario II ---------------------------------------------------------
# 
# 
# ggplot(data = d_rirs$data_short %>% slice_head(n = n),
#        aes(x = id, y = survtime, colour = as.factor(event))) +
#   geom_point()
# ggplot(d_rirs$data, aes(x = obstime, y = y, color = id)) +
#   geom_line() +
#   facet_grid(~marker) +
#   theme(legend.position = "none")
# ggplot(d_rirs$data_hypo %>% filter(event == 1), aes(x = obstime, y = mu, color = id)) +
#   geom_line() +
#   facet_grid(~marker, scales = "free") +
#   theme(legend.position = "none")
# 
# 
surv <- list()
for(i in 100:199) {
  set.seed(i)
  aha <- simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                    probmiss = 0.75, maxfac = 1.75,
                    nmark = 6, long_assoc = "param", M = NULL, 
                    FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                    re_cov_mat = cov,
                    ncovar = 2,
                    lambda = function(t, x) {
                      1.37 * t^(0.37)
                    },
                    gamma = function(x) {
                      -1.5 + 0.48*x[, 3]
                    },
                    alpha = list(function(t, x) {
                      1.5 + 0*t
                    }, function(t, x) {
                      0.6 + 0*t
                    }, function(t, x) {
                      0.3 + 0*t
                    }, function(t, x) {
                      -0.3 + 0*t
                    }, function(t, x) {
                      -0.6 + 0*t
                    }, function(t, x) {
                      -1.5 + 0*t
                    }),
                    mu = list(function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 1] + r[, 2]*t
                    }, function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 3] + r[, 4]*t
                    }, function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 5] + r[, 6]*t
                    }, function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 7] + r[, 8]*t
                    }, function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 9] + r[, 10]*t
                    }, function(t, x, r){
                      0 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
                        r[, 11] + r[, 12]*t
                    }),
                    sigma = function(t, x) {
                      log(0.06) + 0*t
                    }, 
                    tmax = NULL, seed = NULL, 
                    full = TRUE, file = NULL)
  surv[[i]] <- aha
}
summary(sapply(surv, function(x) summary(x)["Median"]))
