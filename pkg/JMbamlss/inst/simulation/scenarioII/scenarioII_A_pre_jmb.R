
# Simulation Scenario II - 22/12/09 ---------------------------------------


# JMbayes estimation
# Simulate ten iterations and then try out different number of splines


# Set up R session --------------------------------------------------------


location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if(location == "server_linux") "./simulation/"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
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
library(JMbamlss)


# Setting for the simulation
start <- 10
stop <- 19
setting <- "scen_II_221209/"


# Objects for all clusters ------------------------------------------------

# Number of individuals and other quantities
n <- 300
argvals <- seq(0, 1, by = 0.01)
x <- seq(0, 1, by = 0.1)

# Random covariance matrix
# Set the eigenvalues but the eigenvectors are random
set.seed(1105)
p <- 6
P <- qr.Q(qr(matrix(rnorm(p^2), ncol = p)))
evals <- c(4, 3, 2, 1, 0.5, 0.2)
cov <- crossprod(P, P*(evals))


# Find spline functions
# Marker1
m1sp1 <- splinefun(x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
m1sp2 <- splinefun(x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
m1sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4))
# Marker2
m2sp1 <- splinefun(x, c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
m2sp2 <- splinefun(x, c(0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
m2sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3))

m1 <- funData(argvals = argvals,
              X = matrix(c(m1sp1(argvals), m1sp2(argvals), m1sp3(argvals)), 
                         nrow = 3, byrow = TRUE))
m2 <- funData(argvals = argvals,
              X = matrix(c(m2sp1(argvals), m2sp2(argvals), m2sp3(argvals)), 
                         nrow = 3, byrow = TRUE))

# True multivariate covariance structure
m <- JMbamlss:::MFPCA_cov(cov = cov, basis_funs = list(m1, m2))



# Data simulation function ------------------------------------------------

simulate_data <- Vectorize(function(i) {
  set.seed(i)
  # Simulate the data
  simdat <- JMbamlss:::simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                       max_obs = 15, probmiss = 0.75, maxfac = 3,
                       nmark = 2, long_assoc = "FPC", M = 6, 
                       FPC_bases = m$functions, FPC_evals = m$values, 
                       #mfpc_args = NULL, re_cov_mat = NULL,
                       ncovar = 2,
                       lambda = function(t, x) {
                         1.65 * t^(0.65)
                       },
                       gamma = function(x) {
                         -3 + 0.3*x[, 3]
                       },
                       alpha = list(function(t, x) {
                         1.1 + 0*t
                       }, function(t, x) {
                         1.1 + 0*t
                       }),
                       mu = list(function(t, x, r){
                         0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                       }, function(t, x, r){
                         0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                       }),
                       sigma = function(t, x) {
                         log(0.06) + 0*t
                       }, 
                       tmax = NULL, seed = NULL, 
                       full = TRUE, file = NULL)
  
  
  save(simdat, file = paste0(results_wd, setting, "/data_pre/d", i, ".Rdata"))
  NULL
})



# Estimate the model for different df -------------------------------------

jmb_df_est <- function(i, df = c(3, 4, 5)) {
  set.seed(i)
  
  # Load the data
  load(paste0(results_wd, setting, "data_pre/d", i, ".Rdata"))
  
  # Get the quantile-based knots for comparability
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                        data = simdat$data)[[1]]$knots
  
  # List of formulas with different degrees of freedom
  f_list <- lapply(df, function (d) {
    as.formula(paste("~ bs(obstime, df =", d,
                     ", Boundary.knots = c(0, 1)) | id"))
  })
  
  # Estimate the model using JMbayes
  simdat_jmb <- simdat$data %>% 
    pivot_wider(names_from = marker, values_from = y)
  
  CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                  data = simdat$data_short %>% filter(marker == "m1"))
  
  lapply(f_list, function(f) {
    lm1 <- lme(m1 ~ obstime * x3, 
               random = f, 
               data = simdat_jmb, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm2 <- lme(m2 ~ obstime * x3, 
               random = f, 
               data = simdat_jmb, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    
    jmb <- jm(CoxFit, list(lm1, lm2), time_var = "obstime",
              n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
              cores = 1, n_chains = 1, 
              GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
              save_random_effects = TRUE)$fit_stats
  })
}



# Analysis ----------------------------------------------------------------


data_sim <- simulate_data(start:stop)
mod_stats <- lapply(start:stop, jmb_df_est, df = c(3, 4))
