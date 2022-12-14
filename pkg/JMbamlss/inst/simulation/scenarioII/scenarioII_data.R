
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
#library(JMbamlss)
devtools::load_all()

# Setting for the simulation
start <- 100
stop <- 199
number_cores <- 5
setting <- "scen_II_221209"
dir.create(paste0(results_wd, setting, "/data"), showWarnings = FALSE)
Sys.time()
sessionInfo()


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
m <- MFPCA_cov(cov = cov, basis_funs = list(m1, m2))


# Simulation function -----------------------------------------------------

parallel_data <- function(i) {
  set.seed(i)
  # Simulate the data
  d_rirs <- simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
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
                       tmax = NULL, seed = 1008, 
                       full = TRUE, file = NULL)
  
  
  save(d_rirs, file = paste0(results_wd, setting, "/data/d", i, ".Rdata"))
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
dat_sim <- mclapply(start:stop, parallel_data, mc.cores = number_cores)

