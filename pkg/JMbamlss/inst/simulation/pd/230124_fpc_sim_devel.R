
# Set Up Simulation: Data -------------------------------------------------

# Simulate censored and uncensored data


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
setting <- "scen_II_230117" # 221209: all data sets are the same
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

# Check how to simulate the data that they are comparable to simscenII
# Number of censorings is relatively small for same censoring as in simscenII


# Simulate the data as in simscenII
d_sim <- JMbamlss:::simCensLongi(nsub = n*100, times = seq(0, 1, by = 0.01), 
                                 max_obs = 15, probmiss = 0.75, maxfac = 3,
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
                                 tmax = NULL, seed = 912, 
                                 full = TRUE, file = NULL)

ggplot(d_sim$data_short, aes(x = survtime)) +
  geom_density()+
  ggtitle("Density of Follow-Up Time (Censoring)") +
  theme_bw() +
  labs(colour = "Event", y = "Density", x = "Follow-Up Time")
table((d_sim$data_short %>% mutate(surv1 = survtime == 1))$surv1) /
  n_high
# FALSE   TRUE 
# 0.3368 0.6632 

# Compare to SimScenII data (takes very long to calculate!)
n_high <- 30000
d_rirs <- simMultiJM(nsub = n_high, times = seq(0, 1, by = 0.01), 
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
ggplot(d_rirs$data_short[1: n_high, ], aes(x = survtime)) +
  geom_density()+
  ggtitle("Density of Follow-Up Time (Censoring)") +
  theme_bw() +
  labs(colour = "Event", y = "Density", x = "Follow-Up Time")
table((d_rirs$data_short[1:n_high, ] %>% mutate(surv1 = survtime == 1))$surv1) /
  n_high
# FALSE      TRUE 
# 0.8016333 0.1983667 

# Simulate the data with higher random censoring to match MJM
d_siH <- JMbamlss:::simCensLongi(nsub = n*100, times = seq(0, 1, by = 0.01), 
                                 max_obs = 15, probmiss = 0.75, maxfac = 1.245,
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
                                 tmax = NULL, seed = 1112, 
                                 full = TRUE, file = NULL)
table((d_siH$data_short %>% mutate(surv1 = survtime == 1))$surv1) /
  n_high
# FALSE      TRUE 
# 0.8010333 0.1989667 
ggplot(d_siH$data_short, aes(x = survtime)) +
  geom_density()+
  ggtitle("Density of Follow-Up Time (Censoring)") +
  theme_bw() +
  labs(y = "Density", x = "Follow-Up Time")

comb <- rbind(d_sim$data_short %>% select(survtime),
              d_siH$data_short %>% select(survtime),
              d_rirs$data_short[1: n_high, ] %>% select(survtime)) %>%
  mutate(type = factor(rep(c("Low Cens", "High Cens", "MJM"), each = n_high)))
ggplot(comb, aes(x = survtime, color = type)) +
  geom_density()+
  ggtitle("Density of Follow-Up Time (Censoring)") +
  theme_bw() +
  labs(colour = "Setting", y = "Density", x = "Follow-Up Time")


# Calculate for one simulation run ----------------------------------------

# Simulate the data as in simscenII
d_sim <- JMbamlss:::simCensLongi(nsub = n, times = seq(0, 1, by = 0.01), 
                                 max_obs = 15, probmiss = 0.75, maxfac = 3,
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
                                 tmax = NULL, seed = 1002, 
                                 full = TRUE, file = NULL)


# Estimate the model using estimated FPCs
# remove observation with less than 3 longitudinal observations
few_obs <- apply(table(d_sim$data_uncens$id, d_sim$data_uncens$marker), 1, 
                 function (x) any(x < 3))
long_obs <- d_sim$data_uncens %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(obstime), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))

mfpca_est <- JMbamlss:::preproc_MFPCA(d_sim$data_uncens %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                                      npc = 3, nbasis = 10,
                                      fit = TRUE, save_uniGAM = TRUE)

# To compare estimated FPCs: Flip, then norm of difference
flipped <- flipFuns(m$functions, mfpca_est$functions)
norm(m$functions - flipped)
sq_diff <- ((m$functions - flipped)^2)
rowSums(sq_diff@.Data[[1]]@X) + rowSums(sq_diff@.Data[[2]]@X)
norm(sq_diff)

# To compare the estimated eigenspace use subspace function
# Conversion from radians to degree: *180/pi
espace_tru <- t(cbind(m$functions$m1@X, m$functions$m2@X))
espace_est <- t(cbind(mfpca_est$functions$m1@X, mfpca_est$functions$m2@X))
espace_fli <- t(cbind(flipped$m1@X, flipped$m2@X))
espace_var <- varimax(espace_fli)$loadings
pracma::subspace(espace_tru, espace_est) * 180/pi
pracma::subspace(espace_fli, espace_est) * 180/pi
pracma::subspace(espace_tru, espace_var) * 180/pi

# To compare the reconstructed observations
uniGAM <- attr(mfpca_est, "uniGAM")
uni_pred <- multiFunData(mapply(function(m, d) {
  p <- predict(m, newdata = d)
  funData(argvals = list(d$obstime[1:101]),
          X = t(matrix(p, nrow = 101, ncol = n)))
}, m = uniGAM, d = split(d_sim$data_full, d_sim$data_full$marker)))
pred <- uni_pred + mfpca_est$fit
tru <- multiFunData(lapply(split(d_sim$data_full, d_sim$data_full$marker), 
                           function (d) {
                             funData(argvals = list(d$obstime[1:101]),
                                     X = t(matrix(d$mu, nrow = 101, ncol = n)))
                           }))
norm(tru-pred)
