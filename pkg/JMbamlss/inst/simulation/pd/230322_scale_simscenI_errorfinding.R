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



# Covariance for Simulation -----------------------------------------------


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
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)),
                                         auto)

# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
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

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Original Data Simulation",
          "Simulated Data Set 1")





# -------------------------------------------------------------------------
# Increase the RE Variance ------------------------------------------------
# -------------------------------------------------------------------------

cov <- cov*20


# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                     max_obs = 15, probmiss = 0.75, maxfac = 1.75,
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

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Scaled Cov, Orig Error Var",
          "Simulated Data Set 1")


# Estimate the model using JMbayes
kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                      data = d_rirs$data)[[1]]$knots
d_rirs_jmb <- d_rirs$data %>% 
  pivot_wider(names_from = marker, values_from = y)
CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                data = d_rirs$data_short %>% filter(marker == "m1"))
lm1 <- lme(m1 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm2 <- lme(m2 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm3 <- lme(m3 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm4 <- lme(m4 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm5 <- lme(m5 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm6 <- lme(m6 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
jmb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "obstime",
          n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
          cores = 1, n_chains = 1, 
          GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
          save_random_effects = TRUE)


# Estimate the model using JMbamlss
seq1 <- seq(0, 1, by = 0.01)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)
mfpca_tru <- JMbamlss:::MFPCA_cov(cov = cov, basis_funs = b_funs)
nfpc <- 12
mfpca_tru_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_tru_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
set.seed(1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, n.iter = 5500,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1)


# -------------------------------------------------------------------------
# Additionally Increase Error Variance ------------------------------------
# -------------------------------------------------------------------------


# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
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
                                  log(0.4) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Scaled Cov, Increased Error Variance",
          "Simulated Data Set 1")

# Estimate the model using JMbamlss
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
set.seed(1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, n.iter =500,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1)



# Increase Error Variance to Mauff Maximum --------------------------------


# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
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
                                  log(0.94) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Scaled Cov, Mauff Error Variance",
          "Simulated Data Set 1")

# Estimate the model using JMbayes
kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                      data = d_rirs$data)[[1]]$knots
d_rirs_jmb <- d_rirs$data %>% 
  pivot_wider(names_from = marker, values_from = y)
CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                data = d_rirs$data_short %>% filter(marker == "m1"))
lm1 <- lme(m1 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm2 <- lme(m2 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm3 <- lme(m3 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm4 <- lme(m4 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm5 <- lme(m5 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm6 <- lme(m6 ~ obstime * x3, random = ~ obstime | id, 
           data = d_rirs_jmb, na.action = na.omit,
           control = lmeControl(opt = "optim"))
jmb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "obstime",
          n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
          cores = 1, n_chains = 1, 
          GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
          save_random_effects = TRUE)

# Estimate the model using JMbamlss
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
set.seed(1)
# Fails after iteration 100
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, n.iter = 5500,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 101, sampler = FALSE,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1, 
                par_trace = TRUE)

alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:101) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 82, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters for Scaled Cov + Mauff Error Var")


# Estimate the model using univariate JM - bamlss
f_uni <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~  obstime:x3 + s(id, bs = "re") + s(obstime, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_uni1 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m1"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE) 
b_uni2 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m2"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni3 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m3"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni4 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m4"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni5 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m5"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni6 <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                   filter(marker == "m6"),
                 timevar = "obstime", idvar = "id", maxit = 5000, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
save(list = paste0("b_uni", 1:6), file = paste0(server_wd, "scen_mauff", 
                                                "/uni/scaleuni_", 1, ".Rdata"))
# all models do not converge
# Only model 5 shows zero acceptance probability for sigma predictor
summary(b_uni5$samples[[1]][, grep("\\.alpha", colnames(b_uni5$samples[[1]]))])
summary(b_uni6$samples[[1]][, grep("\\.alpha", colnames(b_uni6$samples[[1]]))])

table(d_rirs$data$marker, d_rirs$data$id)
b_uni4$parameters$sigma

# Re-start the sampler to see how the sigma traceplots behave
b_uni1short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m1"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni1),
                      optimizer = FALSE, n.iter = 1500,
                      burnin = 0, thin = 1, verbose = TRUE)
plot(b_uni1short, model = "sigma", which = "samples", ask = FALSE)
b_uni2short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m2"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni2),
                      optimizer = FALSE, n.iter = 1500,
                      burnin = 0, thin = 1, verbose = TRUE)
plot(b_uni2short, model = "sigma", which = "samples", ask = FALSE)
b_uni3short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m3"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni3),
                      optimizer = FALSE, n.iter = 1500,
                      burnin = 0, thin = 1, verbose = TRUE)
plot(b_uni3short, model = "sigma", which = "samples", ask = FALSE)
b_uni4short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m4"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni4),
                      optimizer = FALSE, n.iter = 1500,
                      burnin = 0, thin = 1, verbose = TRUE)
plot(b_uni4short, model = "sigma", which = "samples", ask = FALSE)
b_uni5short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m5"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni5),
                      optimizer = FALSE, n.iter = 1500,
                      burnin = 0, thin = 1, verbose = TRUE)
summary(b_uni5short$samples[[1]][, grep("\\.alpha", 
                                        colnames(b_uni5short$samples[[1]]))])
plot(b_uni5short, model = "sigma", which = "samples", ask = FALSE)
b_uni6short <- bamlss(f_uni, family = "jm", data = d_rirs$data %>%
                        filter(marker == "m6"),
                      timevar = "obstime", idvar = "id",
                      start = parameters(b_uni6),
                      optimizer = FALSE, n.iter = 5500,
                      burnin = 0, thin = 1, verbose = TRUE)
summary(b_uni6short$samples[[1]][, grep("\\.alpha", 
                                        colnames(b_uni6short$samples[[1]]))])
plot(b_uni6, model = "sigma", which = "samples", ask = FALSE)


# Re-estimate the model with updating edfs
set.seed(1)
b_uni1_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m1"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE) 
b_uni2_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m2"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE)
b_uni3_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m3"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE)
b_uni4_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m4"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE) 
b_uni5_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m5"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE) 
b_uni6_nu <- bamlss(f_uni, family = "jm", data = d_rirs$data %>% 
                      filter(marker == "m6"),
                    timevar = "obstime", idvar = "id", maxit = 1500,
                    n.iter = 5500, update.nu = TRUE,
                    burnin = 500, thin = 5, verbose = TRUE) 


# Original Covariance and Increased Error Variance ------------------------

cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)),
                                         auto)

# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
                                nmark = 6, long_assoc = "param", M = NULL, 
                                FPC_bases = NULL, FPC_evals = NULL, 
                                mfpc_args = NULL,
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
                                  log(0.94) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Orig Cov, Mauff Error Var - Problem",
          "Simulated Data Set 1")

# Estimate the model using JMbamlss
nfpc <- 6
mfpca_tru <- JMbamlss:::MFPCA_cov(cov = cov, basis_funs = b_funs)
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
f_tru <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_tru_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
mfpca_tru_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
set.seed(1)
# Fails after iteration 236
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, n.iter = 5500,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 237, sampler = FALSE,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1, 
                par_trace = TRUE)

alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:237) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 32, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters for Original Cov + Mauff Error Var")


# Now update nu
set.seed(1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, verbose = TRUE,
                update_nu = TRUE)
summary(b_est$samples[[1]][, grep("accepted", colnames(b_est$samples[[1]]))])
# works

# Increase number of observations
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = 400, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
                                nmark = 6, long_assoc = "param", M = NULL, 
                                FPC_bases = NULL, FPC_evals = NULL, 
                                mfpc_args = NULL,
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
                                  log(0.94) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
set.seed(1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, verbose = TRUE,
                update_nu = TRUE)


# Increase Error Variance to 0.4 ------------------------------------------

# Simulate the data
set.seed(1)
d_rirs <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                                max_obs = 15, probmiss = 0.75, maxfac = 1.75,
                                nmark = 6, long_assoc = "param", M = NULL, 
                                FPC_bases = NULL, FPC_evals = NULL, 
                                mfpc_args = NULL,
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
                                  log(0.4) + 0*t
                                }, 
                                tmax = NULL, seed = NULL, 
                                full = TRUE, file = NULL)

ggplot(d_rirs$data,
       aes(x = obstime, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories - Orig Cov, Incr Error Var",
          "Simulated Data Set 1")

# Estimate the model using JMbamlss
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)
set.seed(1)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, n.iter =500,
                burnin = 500, thin = 5, verbose = TRUE, nu = 1)

# Note: if sigma = function(x) log(0.4) then no problem
