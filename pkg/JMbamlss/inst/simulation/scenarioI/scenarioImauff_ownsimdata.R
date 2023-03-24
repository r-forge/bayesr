
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
server_wd <- switch(location,
                    "workstation" = paste0("/run/user/1000/gvfs/smb-share:",
                                           "server=clapton.wiwi.hu-berlin.de,",
                                           "share=volkmana.hub/JMbamlss/",
                                           "simulation/"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation/",
                    "server_windows" = "H:/JMbamlss/simulation/")


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
library(JMbamlss)


# Setting for the simulation
start <- 100
stop <- 119
number_cores <- 2
setting_old <- "scen_mauff"
setting <- "scen_mauff_2"
# dir.create(paste0(server_wd, setting), showWarnings = FALSE)
# dir.create(paste0(server_wd, setting, "/data"), showWarnings = FALSE)
# dir.create(paste0(server_wd, setting, "/bamlss_tru"), showWarnings = FALSE)
# dir.create(paste0(server_wd, setting, "/bamlss_est"), showWarnings = FALSE)
# dir.create(paste0(server_wd, setting, "/jmb"), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Objects for all clusters ------------------------------------------------

i <- 1
d_mauff <- readRDS(paste0(server_wd, setting_old, "/data/data", i, ".rds"))

# Number of individuals
n <- 150

# Covariance matrix for the data generation

D_num <- matrix(NA, nrow = 12, ncol = 12)
D_num[1:4, 1:4] <- matrix(c(0.919,  0.003,  0.220, -0.008,
                            0.003,  0.006,  0.005, -0.002,
                            0.220,  0.005,  0.421, -0.013,
                            -0.008, -0.002, -0.013,  0.008), byrow = TRUE)
D_num[5:8, 5:8] <- matrix(c(0.551,  0.007, -0.141,  0.015,
                            0.007,  0.014, -0.005, -0.005,
                            -0.141, -0.005,  0.204, -0.042,
                            0.015, -0.005, -0.042,  0.043), byrow = TRUE)  
D_num[9:12, 9:12] <- matrix(c(0.110, -0.028,  0.176,  0.021,
                              -0.028,  0.035, -0.026, -0.003,
                              0.176, -0.026,  3.580,  0.040,
                              0.021, -0.003,  0.040,  0.197), byrow = TRUE)  
D_num[1:4, 5:8] <- matrix(c(0.308,  0.005,  0.265, -0.013,
                            0.012,  0.001,  0.012,  0.002,
                            -0.073, -0.007, -0.172,  0.018,
                            0.007,  0.005,  0.012, -0.013), byrow = TRUE)
D_num[1:4, 9:12] <- matrix(c(0.049,  0.005,  0.124, -0.013,
                             -0.007, -0.005, -0.010,  0.012,
                             0.698,  0.006,  0.680, -0.027,
                             0.056,  0.006,  0.034,  0.004), byrow = TRUE)  
D_num[5:8, 9:12] <- matrix(c(0.095,  0.004, -0.144,  0.032,
                             -0.013,  0.005,  0.035, -0.037,
                             0.826,  0.018, -0.263,  0.024,
                             0.077,  0.019, -0.015,  0.006), byrow = TRUE)
D_num[5:8, 1:4] <- t(D_num[1:4, 5:8])
D_num[9:12, 1:4] <- t(D_num[1:4, 9:12])
D_num[9:12, 5:8] <- t(D_num[5:8, 9:12])

D_numsym <- forceSymmetric(D_num) 
D <- D_numsym



# Estimate true FPCs ------------------------------------------------------


# Use true dependence structure to get a good estimate of the underlying FPC
# basis
# Basis functions on each dimension
seq1 <- seq(0, max(d_mauff$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs)
#save(mfpca_tru, file = paste0(server_wd, setting, "/mfpca_tru.Rdata"))





# Simulate Data -----------------------------------------------------------


set.seed(1)
# Simulate the data
d_rirs <- JMbamlss:::simMultiJM(nsub = n, 
                                times = seq(0, max(d_mauff$Time),
                                            length.out = 101), 
                                probmiss = 0.75, maxfac = 1.75,
                                nmark = 6, long_assoc = "param", M = NULL, 
                                FPC_bases = NULL, FPC_evals = NULL, 
                                mfpc_args = NULL,
                                re_cov_mat = as.matrix(D),
                                ncovar = 2,
                                lambda = function(t, x) {
                                  log(1.65) + 0.65*t
                                },
                                gamma = function(x) {
                                  -5.8 + 0.5*x[, 3]
                                },
                                alpha = list(function(t, x) {
                                  0.1 + 0*t
                                }, function(t, x) {
                                  -0.6 + 0*t
                                }, function(t, x) {
                                  -0.1 + 0*t
                                }, function(t, x) {
                                  -1.41 + 0*t
                                }, function(t, x) {
                                  -1.81 + 0*t
                                }, function(t, x) {
                                  0.75 + 0*t
                                }),
                                mu = list(function(t, x, r){
                                  4.93 + 0.4*t +
                                    r[, 1] + r[, 2]*t
                                }, function(t, x, r){
                                  3.58 + 1*t +
                                    r[, 3] + r[, 4]*t
                                }, function(t, x, r){
                                  1.46 - 0.7*t +
                                    r[, 5] + r[, 6]*t
                                }, function(t, x, r){
                                  1.78 + 0.5*t +
                                    r[, 7] + r[, 8]*t
                                }, function(t, x, r){
                                  0.31 - 1.2*t +
                                    r[, 9] + r[, 10]*t
                                }, function(t, x, r){
                                  1.71 - 0.3*t +
                                    r[, 11] + r[, 12]*t
                                }),
                                sigma = function(t, x) {
                                  log(
                                    I(x$marker == "m1")*0.86 +
                                    I(x$marker == "m2")*0.39 +
                                    I(x$marker == "m3")*0.94 +
                                    I(x$marker == "m4")*0.40 +
                                    I(x$marker == "m5")*0.36 +
                                    I(x$marker == "m6")*0.57
                                  ) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)

save(d_rirs, file = paste0(server_wd, setting, "/data/d", i, ".Rdata"))

# Plot the data
ggplot(data = d_rirs$data_short %>% slice_head(n = n),
       aes(x = id, y = survtime, colour = as.factor(event))) +
  geom_point()
ggplot(d_rirs$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
ggplot(d_rirs$data_hypo %>% filter(event == 1), aes(x = obstime, y = mu, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
summary(c(table(d_rirs$data$id)))
table(d_rirs$data_short[seq_len(n), "event"])


# data with more observations
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = 500, 
                                times = seq(0, max(d_mauff$Time),
                                            length.out = 101), 
                                probmiss = 0.75, maxfac = 1.75,
                                nmark = 6, long_assoc = "param", M = NULL, 
                                FPC_bases = NULL, FPC_evals = NULL, 
                                mfpc_args = NULL,
                                re_cov_mat = as.matrix(D),
                                ncovar = 2,
                                lambda = function(t, x) {
                                  log(1.65) + 0.65*t
                                },
                                gamma = function(x) {
                                  -5.8 + 0.5*x[, 3]
                                },
                                alpha = list(function(t, x) {
                                  0.1 + 0*t
                                }, function(t, x) {
                                  -0.6 + 0*t
                                }, function(t, x) {
                                  -0.1 + 0*t
                                }, function(t, x) {
                                  -1.41 + 0*t
                                }, function(t, x) {
                                  -1.81 + 0*t
                                }, function(t, x) {
                                  0.75 + 0*t
                                }),
                                mu = list(function(t, x, r){
                                  4.93 + 0.4*t +
                                    r[, 1] + r[, 2]*t
                                }, function(t, x, r){
                                  3.58 + 1*t +
                                    r[, 3] + r[, 4]*t
                                }, function(t, x, r){
                                  1.46 - 0.7*t +
                                    r[, 5] + r[, 6]*t
                                }, function(t, x, r){
                                  1.78 + 0.5*t +
                                    r[, 7] + r[, 8]*t
                                }, function(t, x, r){
                                  0.31 - 1.2*t +
                                    r[, 9] + r[, 10]*t
                                }, function(t, x, r){
                                  1.71 - 0.3*t +
                                    r[, 11] + r[, 12]*t
                                }),
                                sigma = function(t, x) {
                                  log(
                                    I(x$marker == "m1")*0.86 +
                                    I(x$marker == "m2")*0.39 +
                                    I(x$marker == "m3")*0.94 +
                                    I(x$marker == "m4")*0.40 +
                                    I(x$marker == "m5")*0.36 +
                                    I(x$marker == "m6")*0.57
                                  ) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)
save(d_rirs, file = paste0(server_wd, setting, "/data/d", i, "moreobs", 
                           ".Rdata"))

# Estimate Model with TRUE MFPCs ------------------------------------------


d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = 12)

f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker +",
    paste0(lapply(seq_along(mfpca_tru$values), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })


b_tru <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                    timevar = "obstime", maxit = 1500, n.iter = 5500,
                    burnin = 500, thin = 5, verbose = TRUE)
# Optimizer does not converge, Updates never accept eta_sigma



# Estimate Model with RE --------------------------------------------------

f_rirs <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker +
    s(id, by = marker, bs = "re") +
    s(id, obstime, by = marker, bs = "re"),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

b_rirs <- bamlss(f_rirs, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                 timevar = "obstime", maxit = 1500, n.iter = 100,
                 burnin =0, thin = 5, verbose = TRUE, par_trace = TRUE)



# Standardize the Longitudinal Outcomes First -----------------------------

# Extract measure of univariate variance
take <- which(!apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
      function (x) any(x < 2)))
mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs$data %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      time = "obstime",
                                      uni_mean = "y ~ 1 + obstime",
                                      npc = 2, nbasis = 4,
                                      save_uniFPCA = TRUE)
univars <- data.frame(marker = factor(paste0("m", 1:6)),
                      ev_std = sqrt(sapply(attr(mfpca_est, "uniFPCA"), 
                                           function(x) sum(x$evalues))))

# Standardize the data
d_rirs_std <- d_rirs$data %>%
  left_join(univars, by = "marker") %>%
  mutate(y = y/ev_std)
save(d_rirs, d_rirs_std, 
     file = paste0(server_wd, setting, "/data/d", i, "std", ".Rdata"))
ggplot(d_rirs_std, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")

# Estimate the MFPC basis
mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs_std %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      time = "obstime",
                                      uni_mean = "y ~ 1 + obstime",
                                      npc = 2, nbasis = 4)
mfpca_list <- lapply(seq_along(mfpca_est$values), 
                     function (i, mfpca = mfpca_est) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })
d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs_std, n = 12,
                                     obstime = "obstime")

# Model for different data
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker +",
    paste0(lapply(seq_along(mfpca_tru$values), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                timevar = "obstime", maxit = 1500, n.iter = 900,
                burnin = 1000, thin = 3, verbose = TRUE)

b_rirs <- bamlss(f_rirs, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                 timevar = "obstime", maxit = 1500, n.iter = 100,
                 burnin =0, thin = 5, verbose = TRUE, par_trace = TRUE)
