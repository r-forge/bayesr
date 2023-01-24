
# Crossvalidation ---------------------------------------------------------


# Set Up ------------------------------------------------------------------


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


setting <- "scen_II_221209/"
dat_setting <- "scen_II_221209/data_pre/"


# Implement Cross-validation ----------------------------------------------

# Load the data
i <- 10
load(paste0(results_wd, dat_setting, "d", i, ".Rdata"))


K <- 10
folds <- data.frame(id = factor(1:300),
                    fold = sample(rep(seq_len(K), 30), size = 300))
d_cv <- left_join(simdat$data, folds, by = "id")

for (k %in% seq_len(K)) {
  
  # Exclude data in the fold
  d_train <- d_cv %>% filter(fold != k)
  d_test <- d_cv %>% filter(fold == k) %>% select(!matches("fpc"))
  
  # Exclude data on other criteria (too few observations, does not span the 
  # time interval)
  few_obs <- apply(table(d_train$id, d_train$marker), 1, function (x) any(x < 3))
  long_obs <- d_train %>%
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
  
  # Calculate MFPCA
  mfpca_est <- JMbamlss:::preproc_MFPCA(
    d_train %>% filter(id %in% take) %>% droplevels(),
    uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
    npc = 3, nbasis = 10)
  
  # Attach FPCs to test data set and account for fixed effects
  d_test <- JMbamlss:::attach_wfpc(mfpca_est, d_test)
  d_test$res <- do.call(c, mapply(function(mod, dat) {
    fit <- predict(mod, newdata = dat, type = "response")
    dat$y - fit
  }, mod = attr(mfpca_est, "uniGAM"), dat = split(d_test, d_test$marker),
  SIMPLIFY = FALSE))
  
}




aha <- as.data.frame(predict(mfpca_est))
oho <- d_train %>% filter(id %in% take) %>% droplevels() %>% 
  select(id, marker, obstime, y) %>%
  mutate(obs = id, argvals1 = obstime, id = NULL, obstime = NULL)
oho <- aha[[1]] %>% left_join(oho, by = c("obs", "obstime"))















# Recreate true MFPCA -----------------------------------------------------

# Number of individuals and other quantities
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



