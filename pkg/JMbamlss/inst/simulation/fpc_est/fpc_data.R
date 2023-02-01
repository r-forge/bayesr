
# Set Up Simulation: Data -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/fpc_estimation"
        else "H:/fpc_estimation")
  results_wd <- if(location == "server_linux") "./simulation/"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/fpc_estimation/simulatio",
                       "n/")
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
setting <- "230124_parallel_to_simscenII"
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
m <- JMbamlss:::MFPCA_cov(cov = cov, basis_funs = list(m1, m2))
m$functions <- multiFunData("m1" = extractObs(m$functions@.Data[[1]]),
                            "m2" = extractObs(m$functions@.Data[[2]]))


# Simulation function -----------------------------------------------------

parallel_data <- function(i) {
  set.seed(i)
  # Simulate the data with higher random censoring to match MJM
  d_sim <- JMbamlss:::simCensLongi(nsub = n, times = seq(0, 1, by = 0.01), 
                                   max_obs = 15, probmiss = 0.75,
                                   maxfac = 1.245,
                                   nmark = 2, long_assoc = "FPC", M = 6, 
                                   FPC_bases = m$functions,
                                   FPC_evals = m$values, ncovar = 2,
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
  
  
  
  saveRDS(d_sim, file = paste0(results_wd, setting, "/data/d", i, ".rds"))
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
dat_sim <- mclapply(start:stop, parallel_data, mc.cores = number_cores)

