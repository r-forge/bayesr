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


# Calculate the true MFPCs ------------------------------------------------

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
m$functions <- multiFunData("m1" = extractObs(m$functions@.Data[[1]]),
                            "m2" = extractObs(m$functions@.Data[[2]]))



# Calculate the MFPCs on each simulation run ------------------

data_fpc_i <- function(i, dat_setting = "scen_II_230117", npc = 3,
                       nbasis = 10, rds = TRUE) {
  
  # Load the data
  if (rds) {
    d_rirs <- readRDS(paste0(results_wd, dat_setting, "/data/d", i, ".rds"))
  } else {
    load(paste0(results_wd, dat_setting, "/data/d", i, ".Rdata"))
  }

  
  # Estimate the model using estimated FPCs
  # remove observation with less than npc longitudinal observations
  few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                   function (x) any(x < npc))
  long_obs <- d_rirs$data %>%
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
  
  mfpca_est <- JMbamlss:::preproc_MFPCA(
    d_rirs$data %>% filter(id %in% take) %>% droplevels(), 
    uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
    npc = npc, nbasis = nbasis)
  
  mfpca_est
}
data_fpc <- Vectorize(data_fpc_i, vectorize.args = "i", SIMPLIFY = FALSE)

# Estimate the FPCs on all simulation settings
mfpc_est4 <- data_fpc(100:199, npc = 4)
mfpc_est3 <- data_fpc(100:199, npc = 3)
mfpc_old <- data_fpc_i(100, dat_setting = "scen_II_221209", npc = 3, 
                       rds = FALSE)
mfpc_est <- c(mfpc_est3, list(mfpc_old))



# Extract Eigenfunctions and Combine to Data.frame ------------------------

eigenfun_data_i <- function(mfpca, true_mfpca) {
  
  flipped <- flipFuns(true_mfpca$functions, mfpca$functions)
  dat <- lapply(flipped@.Data, function(m) {
    data.frame(t = m@argvals[[1]], 
               X = c(t(m@X)), 
               fpc = factor(paste0("fpc.", rep(seq_along(mfpca$values), 
                                               each = length(m@argvals[[1]]))),
                            levels = paste0("fpc.", seq_along(mfpca$values))))
  })
  dat <- do.call(rbind, Map(cbind, 
                            marker = factor(paste0("m", seq_along(dat))), dat))
  dif <- true_mfpca$functions - flipped
  dif_dat <- lapply(dif@.Data, function(m){
    data.frame(t = m@argvals[[1]],
               X = c(t(m@X)),
               fpc = factor(paste0("fpc.", rep(seq_along(mfpca$values), 
                                               each = length(m@argvals[[1]]))),
                            levels = paste0("fpc.", seq_along(mfpca$values))))
  })
  dif_dat <- do.call(rbind, Map(cbind, 
                                marker = factor(paste0("m", 
                                                       seq_along(dif_dat))), 
                                dif_dat))
  d <- rbind(dat, dif_dat)
  d$type <- factor(rep(c("abs", "dif"), c(nrow(dat), nrow(dif_dat))))
  d
  
}

eigenfun_data_list <- Vectorize(eigenfun_data_i, vectorize.args = "mfpca",
                                SIMPLIFY = FALSE)

eigenfun_data <- function(mfpca_list) {
  dat <- eigenfun_data_list(mfpca_list, m)
  dat <- do.call(rbind, Map(cbind, it = factor(seq_along(dat)), dat))
  dat_true <- lapply(m$functions@.Data, function(mt){
    data.frame(t = mt@argvals[[1]],
               X = c(t(mt@X)),
               fpc = factor(paste0("fpc.", rep(seq_along(m$values), 
                                               each = length(mt@argvals[[1]]))),
                            levels = paste0("fpc.", seq_along(m$values))),
               type = factor("abs", levels = c("abs", "dif")),
               it = factor(0))
  })
  dat_true <- do.call(rbind, Map(cbind, 
                                 marker = factor(paste0("m", 
                                                        seq_along(dat_true))), 
                                 dat_true))
  rbind(dat, dat_true)
  
}

dat <- eigenfun_data(mfpc_est)



# Find best fitting FPC estimate ------------------------------------------

# Over all FPCs 
# Best estimates
d <- dat %>%
  filter(type == "dif") %>%
  group_by(it) %>%
  summarize(SE = sum(X^2)) %>%
  arrange(SE)  %>%
  mutate(type = factor("10"))


# Worst estimates
dat %>%
  filter(type == "dif") %>%
  group_by(it) %>%
  summarize(SE = sum(X^2)) %>%
  arrange(SE) %>%
  tail(n = 10)

# Iteration 1 among the top 10 worst estimates

# Over the first three FPCs
dat %>%
  filter(type == "dif", fpc %in% c("fpc.1", "fpc.2", "fpc.3")) %>%
  group_by(it) %>%
  summarize(SE = sum(X^2)) %>%
  arrange(SE) %>%
  head(n = 10)

# Iteration 64 is second best estimate based on squared error of estimated FPCs
# over all FPCs and third best over the first three


# Plot estimated vs True FPC ----------------------------------------------

ggplot(dat %>% 
         filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3"), type == "abs") %>%
         mutate(orig = factor(case_when(it == 0 ~ "True",
                                        it == 101 ~ "It 1008", 
                                        it == 64 ~ "It 64",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("True", "It 1008", "It 64", "Other")),
                seq = case_when(it == 0 ~ 105, it == 101 ~ 102, it == 64 ~ 103,
                                TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig, 
           size = orig)) +
  geom_line() +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 1, 0.2)) +
  scale_color_manual(values = c("red", "orange", "cornflowerblue", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Estimated Eigenfunctions and True Eigenfunction")

ggplot(dat %>% 
         filter(fpc %in% c("fpc.4", "fpc.5", "fpc.6"), type == "abs") %>%
         mutate(orig = factor(case_when(it == 0 ~ "True",
                                        it == 101 ~ "It 1008", 
                                        it == 64 ~ "It 64",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("True", "It 1008", "It 64", "Other")),
                seq = case_when(it == 0 ~ 105, it == 101 ~ 102, it == 64 ~ 103,
                                TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig, 
           size = orig)) +
  geom_line() +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 1, 0.2)) +
  scale_color_manual(values = c("red", "orange", "cornflowerblue", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Estimated Eigenfunctions and True Eigenfunction")



# Plot differences between estimate and true FPC --------------------------

ggplot(dat %>% 
         filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3"), type == "dif") %>%
         mutate(orig = factor(case_when(it == 101 ~ "It 1008", 
                                        it == 64 ~ "It 64",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("It 1008", "It 64", "Other")),
                seq = case_when(it == 0 ~ 105, it == 101 ~ 102, it == 64 ~ 103,
                                TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig,
           size = orig)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 0.2)) +
  scale_color_manual(values = c("orange", "cornflowerblue", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Difference between True and Estimated FPC")

ggplot(dat %>% 
         filter(fpc %in% c("fpc.4", "fpc.5", "fpc.6"), type == "dif") %>%
         mutate(orig = factor(case_when(it == 101 ~ "It 1008", 
                                        it == 64 ~ "It 64",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("It 1008", "It 64", "Other")),
                seq = case_when(it == 0 ~ 105, it == 101 ~ 102, it == 64 ~ 103,
                                TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig, 
           size = orig)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 0.2)) +
  scale_color_manual(values = c("orange", "cornflowerblue", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Difference between True and Estimated FPC")



# Estimation of FPCs with Different Basis Functions -------------------------


# Estimate the FPCs using only four basis functions
mfpc_est_4b <- data_fpc(100:199, npc = 3, nbasis = 4)
mfpc_old_4b <- data_fpc_i(100, dat_setting = "scen_II_221209", npc = 3,
                          nbasis = 4, rds = FALSE)
mfpc_est_4b <- c(mfpc_est_4b, list(mfpc_old_4b))
dat_4b <- eigenfun_data(mfpc_est_4b)

# Best estimates
d_4b <- dat_4b %>%
  filter(type == "dif") %>%
  group_by(it) %>%
  summarize(SE = sum(X^2)) %>%
  arrange(SE) %>%
  mutate(type = factor("4"))

# Estimate the FPCs using only twenty basis functions
mfpc_est_20b <- data_fpc(100:199, npc = 3, nbasis = 20)
mfpc_old_20b <- data_fpc_i(100, dat_setting = "scen_II_221209", npc = 3,
                          nbasis = 20, rds = FALSE)
mfpc_est_20b <- c(mfpc_est_20b, list(mfpc_old_20b))
dat_20b <- eigenfun_data(mfpc_est_20b)

# Best estimates
d_20b <- dat_20b %>%
  filter(type == "dif") %>%
  group_by(it) %>%
  summarize(SE = sum(X^2)) %>%
  arrange(SE) %>%
  mutate(type = factor("20"))

# Compare the Sum of Squares
dat_comb <- rbind(d, d_4b, d_20b) %>%
  mutate(type = factor(type, levels = c("4", "10", "20")))
ggplot(dat_comb, aes(x = type, y = SE)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Number of Basis Functions", y = "SSE") +
  ggtitle("Sum of Squared Errors in Estimated MFPCs",
          paste("Different Number of Basis Functions in Estimation of", 
                "Univariate Covariance"))

# Compare the estimate of 4 basis functions with truth
ggplot(dat_4b %>% 
         filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3"), type == "abs") %>%
         mutate(orig = factor(case_when(it == 0 ~ "True",
                                        it == 101 ~ "It 1008",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("True", "It 1008", "Other")),
                seq = case_when(it == 0 ~ 105, TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig, 
           size = orig)) +
  geom_line() +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 0.2)) +
  scale_color_manual(values = c("red", "orange", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Estimated Eigenfunctions and True Eigenfunction",
          "Univariate FPCA with 4 Basis Functions")

# Compare the estimate of 20 basis functions with truth
ggplot(dat_20b %>% 
         filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3"), type == "abs") %>%
         mutate(orig = factor(case_when(it == 0 ~ "True",
                                        it == 101 ~ "It 1008",
                                        it %in% as.character(1:100) ~ "Other"),
                              levels = c("True", "It 1008", "Other")),
                seq = case_when(it == 0 ~ 105, TRUE ~ as.numeric(it))),
       aes(x = t, y = X, group = seq, color = orig, alpha = orig, 
           size = orig)) +
  geom_line() +
  theme_bw() +
  scale_alpha_manual(values = c(1, 1, 0.2)) +
  scale_color_manual(values = c("red", "orange", "black")) +
  scale_size_manual(values = c(1, 0.8, 0.5)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = paste(c(0, 0.25, 0.5, 0.75, 1))) +
  facet_grid(fpc ~ marker) +
  labs(scale = NULL, alpha = NULL, color = NULL, size = NULL) +
  ggtitle("Estimated Eigenfunctions and True Eigenfunction",
          "Univariate FPCA with 20 Basis Functions")



# Angle between the Estimated and True MFPCs ------------------------------


sapply(mfpc_est, function (est, tru = m) {
  ebase <- do.call(rbind,
                   lapply(est$functions@.Data, function (x) t(x@X)))
  tbase <- do.call(rbind,
                   lapply(tru$functions@.Data, function (x) t(x@X)))
  pracma::subspace(ebase, tbase)
})


# Extract NFPC and Ratio of Explained Var Per Dim -------------------------

eigenfun_npc <- Vectorize(function(mfpca) {
  
  vals <- which(mfpca$values > 0)
  nfpc <- min(which(
    cumsum(mfpca$values[vals])/sum(mfpca$values[vals]) > 0.95))
  uni_norms <- lapply(mfpca$functions, norm)
  uni_vars <- lapply(uni_norms, function (m) {
    m*mfpca$values
  })
  tot_vars <- sapply(uni_vars, sum)
  exp_vars <- sapply(uni_vars, function(x) sum(x[seq_len(nfpc)]))
  tot_var <- sum(tot_vars)
  
  data.frame(type = factor(c("nfpc", "tot", rep("tot_vars", 6), 
                             rep("exp_vars", 6), rep("rate_vars", 6))), 
             value = c(nfpc, tot_var, tot_vars, exp_vars, exp_vars/tot_vars), 
             marker = factor(c("all", "all", paste0("m", 1:6), 
                               paste0("m", 1:6), paste0("m", 1:6))))
}, SIMPLIFY = FALSE)

eigenfun_npc_info <- function(mfpca_list) {
  dat <- eigenfun_npc(mfpca_list)
  do.call(rbind, Map(cbind, it = factor(seq_along(dat)), dat))
}

dat <- eigenfun_npc_info(mfpc_est) %>%
  mutate(Problem = it %in% probl)

ggplot(dat %>% filter(type == "nfpc"), aes(x = factor(value), fill = Problem)) + 
  geom_bar(position = position_dodge2(preserve = "single")) +
  theme_bw() +
  ggtitle("Number of FPCs in the Est_95 Models") +
  labs(x = "Number of FPCs", y = "Number of Simulation Iterations")

ggplot(dat %>% filter(type == "tot"), aes(y = log(value), x = Problem, 
                                          fill = Problem)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Log(Total Variance)") +
  ggtitle("Total Variance over Simulations", "Sum of Multivariate Eigenvalues")

ggplot(dat %>% filter(type == "tot_vars"), 
       aes(x = marker, y = log(value), group = it, color = factor(Problem),
           alpha = factor(Problem))) +
  geom_line() +
  theme_bw() +
  labs(y = "Log(Total Univariate Variance)", color = "Problem", 
       alpha = "Problem") +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_color_manual(values = c("black", "orange")) +
  ggtitle("Total Univariate Variance over Simulations", 
          "Sum of Univariate Norms times Multivariate Eigenvalues")

ggplot(dat %>% filter(type == "rate_vars"), 
       aes(x = marker, y = value, group = it, color = factor(Problem),
           alpha = factor(Problem))) +
  geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 1, color = "blue") +
  theme_bw() +
  labs(y = "Explained Univariate Variance", color = "Problem", 
       alpha = "Problem") +
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_color_manual(values = c("black", "orange")) +
  ggtitle("Rate of Explained Univariate Variance over Simulations")



# Compare Level of Explained Variance with MSE ----------------------------
load(paste0(results_wd, "results.Rdata"))

r <- do.call(rbind, Map(cbind, it = factor(seq_along(r_est_95)), r_est_95)) %>%
  filter(type == "MSE", predictor == "mu") %>%
  mutate(marker = factor(marker)) %>%
  group_by(it) %>%
  mutate(max_mse = max(value),
         large_mse = sum(value > 0.1)) %>%
  ungroup()

ggplot(r, aes(x = marker, y = value, group = it, 
              color = factor(max_mse > 0.1), alpha = factor(max_mse > 0.1))) +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.1, 0.5)) +
  labs(y = "MSE", alpha = "Max(MSE)>0.1", color = "Max(MSE)>0.1") +
  ggtitle("MSE of mu for Iterations with Maximal MSE > 0.1")
  
ggplot(r %>% group_by(it) %>% slice_head(n = 1), aes(x = large_mse)) +
  geom_bar() +
  theme_bw() +
  labs(x = "Number of MSE > 0.1", y = "Number of Simulation Iterations") +
  ggtitle("Number of Large MSE of mu Values")


d <- dat %>% 
  filter(type == "rate_vars") %>%
  left_join(r, by = c("it", "marker"))

ggplot(d, aes(x = value.x, y = value.y)) +
  geom_point() +
  facet_wrap(marker~.) +
  theme_bw() +
  labs(x = "Rate of Explained Variance", y = "MSE") +
  ggtitle("Rate of Explained Variance vs. MSE of mu")


d <- dat %>%
  filter(type == "nfpc") %>%
  right_join(r, by = c("it"))
ggplot(d, aes(x = marker.y, y = value.y, group = it)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~value.x) + 
  labs(y = "MSE", x = "Marker") +
  theme_bw() + 
  ggtitle("Number of FPCs and MSE of mu")




# Priori and Posteriori Eigenvalues ---------------------------------------

# Extract Priori Eigenvalues

extract_eval_info <- Vectorize(function (mfpca) {
  
  # Calculate the number of included eigenvalues
  vals <- which(mfpca$values > 0)
  nfpc <- min(which(
    cumsum(mfpca$values[vals])/sum(mfpca$values[vals]) > 0.95))
  
  data.frame(ev = mfpca$values[seq_len(nfpc)],
             fpc = factor(paste0("fpc.", seq_len(nfpc))))
  
}, SIMPLIFY =  FALSE)

priori_ev <- extract_eval_info(mfpc_est)
pri_ev <- do.call(rbind, Map(cbind, it = seq_along(priori_ev), priori_ev))


# Extract Posteriori 'Eigenvalues'
extract_fpc_info_i <- function (i, setting, folder) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  samples(b_est, model = "mu")
  
}
extract_fpc_info <- Vectorize(extract_fpc_info_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
fpc_info <- function (setting, folder) {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_fpc_info(m, setting = setting, folder = folder)
  
}

b_est_95 <- fpc_info(setting = "scen_I_130922", folder = "bamlss_est")

extract_eval_sample <- Vectorize(function (samp) {

  # Extract the relevant columns
  names <- colnames(samp)
  tau21 <- samp[, grep("\\.tau21", names)]
  names <- names[grep("\\.tau21", names)]
  names <- gsub("mu\\.s\\.s\\(id,", "", names)
  names <- gsub("\\)\\.tau21", "", names)
  tau21 <- matrix(tau21, ncol = length(names))
  colnames(tau21) <- names
  
  data.frame(ev = colMeans(tau21),
             fpc = factor(colnames(tau21)))
}, SIMPLIFY = FALSE)
posteriori_ev <- extract_eval_sample(b_est_95)
pos_ev <- do.call(rbind, Map(cbind, it = seq_along(posteriori_ev), 
                             posteriori_ev))

p <- pri_ev %>%
  filter(! it %in% probl) %>%
  left_join(pos_ev, by = c("it", "fpc"))
ggplot(p, aes(x = log(ev.x), y = log(ev.y))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~fpc) +
  theme_bw() +
  labs(x = "Prior", y = "Posterior") +
  ggtitle("Eigenvalues of Preliminary (Prior) and Posterior",
          "Log(Eigenvalues)")


str(pos_ev)
pos <- pos_ev %>%
  group_by(it) %>%
  arrange(desc(ev)) %>%
  mutate(order = seq_len(n())) %>%
  ungroup() %>%
  arrange(it, fpc)
ggplot(pos, aes(x = order)) +
  geom_bar() +
  facet_wrap(~fpc) +
  theme_bw() +
  labs(y = "Number of Simulation Runs", x = "Order of Eigenvalue") +
  ggtitle("Sequence of Eigenvalues Changes in Posterior")


# Vielleicht die Eigenwerte als Priori-Werte verwenden?



# Correlation of Posterior Sample -----------------------------------------

post_covarianc_mat <- Vectorize(function (samp) {
  
  samp <- samp[, grep("\\.b[0-9]+", colnames(samp))]
  n_fpc <- max(as.numeric(gsub("mu\\.s\\.s\\(id,fpc\\.([0-9]+)\\)\\.b[0-9]+", 
                               "\\1", colnames(samp))))
  cov_mats <- apply(samp, 1, function(row) {
    cov(matrix(row, ncol = n_fpc, nrow = 150))
  })
  rowMeans(matrix(cov_mats, nrow = n_fpc^2))
  
}, SIMPLIFY = FALSE)

post_cov <- post_covarianc_mat(b_est_95)

p_cov <- lapply(post_cov, function (it) {
  
  nfpc <- sqrt(length(it))
  if(nfpc == 1) {
    data.frame(cor = NA, nfpc = nfpc)
  } else {
    data.frame(cor = cov2cor(matrix(it, ncol = nfpc))[
      lower.tri(matrix(it, ncol = nfpc))],
               nfpc = nfpc)
  }
  
})
p_cov <- do.call(rbind, Map(cbind, it = seq_along(p_cov), p_cov))

ggplot(p_cov, aes(y = cor, x = 1)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~nfpc) +
  theme_bw() +
  labs(x = "All Correlation Values over All Simulations", y = "Correlation") +
  ggtitle("Correlation between Random Effects Scatters around 0 ")
