

# Simulation Scenario II - 22/12/09 ---------------------------------------




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

setting <- "230124_parallel_to_simscenII/"
model_type <- c("A_uncens/", "A_cens/", "A_jm/")



# Reconstruct True MFPC ---------------------------------------------------

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


# Functions for Evaluation ------------------------------------------------

angle <- function (A, B) 
{
  if (!is.numeric(A) || !is.numeric(B)) 
    stop("Arguments 'A' and 'B' must be numeric matrices.")
  if (is.vector(A)) 
    A <- matrix(c(A), nrow = length(A), ncol = 1)
  if (is.vector(B)) 
    B <- matrix(c(B), nrow = length(B), ncol = 1)
  if (nrow(A) != nrow(B)) 
    stop("Matrices 'A' and 'B' must have the same number of rows.")
  A <- pracma::orth(A)
  B <- pracma::orth(B)
  if (ncol(A) < ncol(B)) {
    tmp <- A
    A <- B
    B <- tmp
  }
  for (k in 1:ncol(A)) {
    B <- B - A[, k] %*% t(A[, k]) %*% B
  }
  acos(Reduce("*", cos(asin(sapply(svd(B)$d, function (i) min(1, i))))))
}

eval_i <- function(it, mtype, setting, dat_setting, m_true, model) {
  
  mfpc <- readRDS(paste0(results_wd, setting, mtype, "/m", it, ".rds"))
  d_sim <- readRDS(paste0(results_wd, dat_setting, "/data/d", it, ".rds"))
  
  # Which observation were removed in mfpca
  d_sim_rm <- d_sim[[switch(mtype, "uncens" = "data_uncens", 
                            "cens" = "data_cens", "jm" = "data")]]
  few_obs <- apply(table(d_sim_rm$id, d_sim_rm$marker), 1, 
                   function (x) any(x < 3))
  long_obs <- d_sim_rm %>%
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
  
  d_sim <- d_sim[[if (mtype != "jm") "data_full" else "data_hypo"]]
  
  # Evaluate estimated FPCs ---------------------------
  flipped <- flipFuns(m_true, mfpc$functions)
  fpc_norms <- norm(m_true - flipped)
  
  # Evaluate estimated eigenspace ---------------------
  espace_tru <- t(do.call(cbind, lapply(m_true, function(m) m@X)))
  espace_est <- t(do.call(cbind, lapply(mfpc$functions, function(m) m@X)))
  angles <- pracma::subspace(espace_tru, espace_est) * 180/pi
  tot_angle <- angle(espace_tru, espace_est) * 180/pi
  
  # Evaluate reconstructed observations ---------------
  
  # Predict the univariate means
  uni_pred <- multiFunData(mapply(function(m, d) {
    p <- predict(m, newdata = d %>% filter(id %in% take))
    funData(argvals = list(d$obstime[1:101]),
            X = t(matrix(p, nrow = 101, ncol = length(take))))
  }, m = attr(mfpc, "uniGAM"), d = split(d_sim, d_sim$marker)))
  pred <- uni_pred + mfpc$fit
  tru <- multiFunData(
    lapply(split(d_sim, d_sim$marker), 
           function (d) {
             d <- d %>% filter(id %in% take)
             funData(argvals = list(d$obstime[1:101]),
                     X = t(matrix(d$mu, nrow = 101, ncol = length(take))))
             }))
  obs_norms <- norm(tru - pred)
  
  # Combine to data set
  data.frame(it = it,
             type = factor(rep(c("fpc", "angles", "tot_angle", "obs"),
                               times = c(length(fpc_norms), length(angles),
                                         length(tot_angle), 
                                         length(obs_norms)))),
             value = c(fpc_norms, angles, tot_angle, obs_norms),
             seq = c(seq_along(fpc_norms), seq_along(angles),
                     seq_along(tot_angle), seq_along(obs_norms)))
  
}
eval_fpc_sim <- Vectorize(eval_i, vectorize.args = "it", SIMPLIFY = FALSE)


dat_uncens <- do.call(
  rbind,
  eval_fpc_sim(100:199, mtype = "uncens", 
               setting = "230124_parallel_to_simscenII/A_",
               dat_setting = "230124_parallel_to_simscenII/", 
               m_true = m$functions)) %>%
  mutate(setting = factor("uncens"))
dat_cens <- do.call(
  rbind,
  eval_fpc_sim(100:199, mtype = "cens", 
               setting = "230124_parallel_to_simscenII/A_",
               dat_setting = "230124_parallel_to_simscenII/", 
               m_true = m$functions)) %>%
  mutate(setting = factor("cens"))
dat_jm <-do.call(
  rbind,
  eval_fpc_sim(100:199, mtype = "jm", 
               setting = "230124_parallel_to_simscenII/A_",
               dat_setting =  "../../JMbamlss/simulation/scen_II_230117", 
               m_true = m$functions)) %>%
  mutate(setting = factor("jm"))

dat <- rbind(dat_uncens, dat_cens, dat_jm)
saveRDS(dat, file = paste0(results_wd, setting, "eval_dat.rds"))

ggplot(data = dat %>% filter(type == "fpc"), 
       aes(x = factor(seq), y = value, color = setting)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "FPC", y = expression("||" * psi ~ - ~ hat(psi) * "||"^2),
       color = "Setting") +
  ggtitle("Estimation of FPCs", 
          "Norms of Differences Between True and Estimated FPCs")

ggplot(data = dat %>% filter(type == "angles"), 
       aes(y = value, x = setting, color = setting)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = "Angle (in °)", color = "Setting") +
  ggtitle("Angle Between True and Estimated Basis Space")

ggplot(data = dat %>% filter(type == "tot_angle"), 
       aes(y = value, x = setting, color = setting)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = "Angle (in °)", color = "Setting") +
  ggtitle("Totatl Angle Between True and Estimated Basis Space")

ggplot(data = dat %>% filter(type == "obs") %>% group_by(setting, it) %>%
         summarize(Median = median(value)) %>% ungroup(), 
       aes(x = setting, y = Median, color = setting)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = expression("Med("~"||" * X[i] ~ - ~ hat(X[i]) * "||"^2~")"),
       color = "Setting") +
  ggtitle("Median Reconstruction Error of Observations", 
          "Norms of Differences Between True and Fitted Observations")

ggplot(data = dat %>% filter(type == "obs") %>% group_by(setting, it) %>%
         summarize(Sum = sum(value)) %>% ungroup(), 
       aes(x = setting, y = log(Sum), color = setting)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, 
       y = expression("Log(Sum("~"||" * X[i] ~ - ~ hat(X[i]) * "||"^2~"))"),
       color = "Setting") +
  ggtitle("Log(Sum) of All Reconstruction Errors of Observations", 
          "Norms of Differences Between True and Fitted Observations")

ggplot(data = dat %>% filter(type != "obs") %>% 
         pivot_wider(id_cols = c("it", "setting"), names_from = type,
                     values_fn = sum),
       aes(x = fpc, y = angles)) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Sum of FPC Differences", y = "Angles") +
  ggtitle("Agreement Between Angles and Estimation of FPCs",
          paste0("Sum over All Six Norms of Differences of True and Estimated",
                 "FPCs"))


# Compare the different criteria
dat_sum <- dat %>%
  group_by(setting, it, type) %>%
  summarise(Sum = sum(value), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(id_cols = c("setting", "it"), names_from = type, 
              values_from = Sum)

dat_med <- dat %>%
  filter(type == "obs") %>%
  group_by(setting, it) %>%
  summarise(obs_med = median(value), .groups = "keep") %>%
  ungroup() %>%
  right_join(dat_sum, by = c("setting", "it"))


ggplot(data = dat_sum, aes(x = tot_angle, y = angles)) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Total Angle", y = "Angles") +
  ggtitle("Agreement Between Angles and Estimation of FPCs",
          paste0("Sum over All Six Norms of Differences of True and Estimated",
                 "FPCs"))

ggplot(data = dat_sum, aes(x = fpc, y = angles)) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Sum of FPC Differences", y = "Angles") +
  ggtitle("Agreement Between Angles and Estimation of FPCs",
          paste0("Sum over All Six Norms of Differences of True and Estimated",
                 "FPCs"))

ggplot(data = dat_sum, aes(x = fpc, y = log(obs))) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Sum of FPC Differences", y = "Log(Obs)") +
  ggtitle("Agreement Between Reconstructed Observations and Estimation of FPCs",
          "Log-Scale for Sum over Errors For Observations")

ggplot(data = dat_sum, aes(x = angles, y = log(obs))) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Angles", y = "Log(Obs)") +
  ggtitle("Agreement Between Angles and Reconstructed Observations",
          "Log-Scale for Sum over Errors For Observations")
# Sum of FPC Differences seem to agree better

# Compare to median of observations
ggplot(data = dat_med, aes(x = fpc, y = obs_med)) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Sum of FPC Differences", y = "Med(Obs)") +
  ggtitle("Agreement Between Reconstructed Observations and Estimation of FPCs",
          "Median over Errors For Observations")
ggplot(data = dat_med, aes(x = angles, y = obs_med)) +
  geom_point() +
  facet_wrap(~setting, scales = "free") +
  labs(x = "Angles", y = "Med(Obs)") +
  ggtitle("Agreement Between Angles and Reconstructed Observations",
          "Median over Errors For Observations")



# Use MSE from MJM Simulation to Find Best Criterion ----------------------

# Not reliable as F contains incorrect model fits
mjm_fits <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(results_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "F/", data_wd = "data/", name = "JM", rds = TRUE)
saveRDS(mjm_fits, 
        file = paste0(results_wd, "../../JMbamlss/simulation/scen_II_230117/",
                      "simscenII_F.rds"))

# Instead use ...
# readRDS()

mu_mse <- mjm_fits %>% 
  filter(type == "MSE", predictor == "mu") %>%
  group_by(it) %>%
  summarize(Sum = sum(value), .groups = "drop_last") %>%
  arrange(desc(Sum)) %>%
  mutate(it = as.numeric(substr(it, 2, 4)))

final_fit <- dat_med %>%
  filter(setting == "jm") %>%
  mutate(logsumobs = log(obs), obs = NULL) %>%
  pivot_longer(3:6, names_to = "crit") %>%
  left_join(mu_mse, by = "it")

ggplot(final_fit, aes(x = value, y = Sum)) + 
  geom_point() +
  facet_wrap(~crit, scales = "free") +
  theme_bw() +
  labs(x = "Value of Criterion", "Sum of Longitudinal MSEs") +
  ggtitle("Accordance of Different MFPC Criteria with MSE from MJM Fit")
