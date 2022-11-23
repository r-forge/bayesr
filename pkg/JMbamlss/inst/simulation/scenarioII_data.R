
# Set Up Simulation: Data -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
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
start <- 100
stop <- 199
number_cores <- 5
setting <- "scen_I_130922"
dir.create(paste0("simulation/", setting, "/data"), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Objects for all clusters ------------------------------------------------

# Number of individuals
n <- 500

# Covariance matrix for the data generation
cov <- matrix(c(0.92, 0, 0.22, -0.01, 0.31, 0.01, -0.07, 0.01, 0.05, -0.01, 0.7, 0.06,
                0, 0.01, 0, 0, 0, 0, -0.01, 0.01, 0.01, 0, 0.01, 0.01,
                0.22, 0, 0.42, -0.01, 0.26, 0.01, -0.17, 0.01, 0.12, -0.01, 0.68, 0.03,
                -0.01, 0, -0.01,  0.01, -0.01, 0, 0.02, -0.01, -0.01, 0.01, -0.03, 0,
                0.31, 0, 0.26, -0.01, 0.55, 0.01, -0.14, 0.02, 0.1, -0.01, 0.83, 0.08,
                0.01, 0, 0.01, 0, 0.01, 0.01, 0, 0, 0, 0, 0.02, 0.02,
                -0.07, -0.01, -0.17, 0.02, -0.14, 0, 0.2, -0.04, -0.14, 0.03, -0.26, -0.02,
                0.01, 0.01, 0.01, -0.01, 0.02, 0, -0.04, 0.04, 0.03, -0.04, 0.02, 0.01,
                0.05, 0.01, 0.12, -0.01, 0.1, 0, -0.14, 0.03, 0.11, -0.03, 0.18, 0.02,
                -0.01, 0, -0.01, 0.01, -0.01, 0, 0.03, -0.04, -0.03, 0.04, -0.03, 0,
                0.7, 0.01, 0.68, -0.03, 0.83, 0.02, -0.26, 0.02, 0.18, -0.03, 3.58, 0.04, 
                0.06, 0.01, 0.03, 0, 0.08, 0.02, -0.02, 0.01, 0.02, 0, 0.04, 0.2),
              nrow = 12, ncol = 12)


# Simulation function -----------------------------------------------------

parallel_data <- function(i) {
  set.seed(i)
  
  # Simulate the data
  d_rirs <- simMultiJM(nsub = n, times = seq(0, 25, by = 0.2), 
                       max_obs = 15, probmiss = 0.75, maxfac = 1.75,
                       nmark = 6, long_assoc = "param", M = NULL, 
                       FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                       re_cov_mat = cov,
                       ncovar = 2,
                       lambda = function(t, x) {
                         1.65 * t^(0.65)
                       },
                       gamma = function(x) {
                         -1.5 + 0.48*x[, 3]
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
                         4.93 + 0.2*t - 0.25*x[, 3] - 0.05*t*x[, 3] +
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
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_data, mc.cores = number_cores)

