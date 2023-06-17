
# Collinearity in SimScen II ----------------------------------------------

# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}
server_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.", 
                    "hu-berlin.de,share=volkmana.hub")
dat_setting <- "scen_II_230117"


# Always
library(survival)
library(JMbayes2)
library(tidyverse)
devtools::load_all(".")


# Calculate Iteration 100 with Truncation (A2) ----------------------------

# Load the data
load(file.path(server_wd, "JMbamlss/simulation", dat_setting, 
               "data/d100.Rdata"))

# Estimate the model using estimated FPCs
# remove observation with less than 3 longitudinal observations
few_obs <- apply(table(d_sim$data$id, d_sim$data$marker), 1, 
                 function (x) any(x < 3))
long_obs <- d_sim$data %>%
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
  d_sim$data %>% filter(id %in% take) %>% droplevels(), 
  uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
  npc = 3, nbasis = 10)
vals <- which(mfpca_est$values > 0)

# Number of FPCs based on univariate explained variation >95
uni_norms <- lapply(mfpca_est$functions, norm)
nfpc <- max(sapply(uni_norms, function (n) {
  min(which(
    cumsum(mfpca_est$values[vals] * n) / sum(mfpca_est$values[vals] * n) >
      0.95
  ))
}))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i,
                                                        mfpca = mfpca_est){
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
d_sim_est <- d_sim$data[, -grep("fpc\\.", names(d_sim$data))]
d_sim_est <- JMbamlss:::attach_wfpc(mfpca_est, d_sim_est, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps",
                                           xt = list("scale" = FALSE)),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_sim_est,
                timevar = "obstime", maxit = 1500, sampler = FALSE, coll = TRUE,
                verbose = TRUE)


# Analyze Collinearity Model Truncated ------------------------------------

nmarker <- length(mfpca_est$functions)
start <- max(which(sapply(b_est$model.stats$optimizer$coll, is.null))) + 1
stop <- length(b_est$model.stats$optimizer$coll)

# Use data matrix
col_x <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  X <- matrix(x$X, ncol = nmarker)
  svd(crossprod(X))$d
})
condi_x <- apply(col_x, 2, function (x) x[1] / x[nmarker])
ratio_x <- apply(col_x, 2, function (x) x / sum(x))

# Use Hesse matrix
col_i <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  H <- matrix(x$I, ncol = nmarker)
  s <- diag(diag(H)^(-1/2))
  svd(s %*% H %*% s)$d
})
condi_i <- apply(col_i, 2, function (x) x[1] / x[nmarker])
ratio_i <- apply(col_i, 2, function (x) x / sum(x))

condi_dat <- data.frame(
  it = rep(start:stop, times = 2),
  vals = c(condi_x, condi_i),
  type = factor(rep(c(0, 1), each = length(start:stop)), labels = c("XTX", "I_S"))
)
ggplot(condi_dat, aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Condition Numbers for SimScen II It 100 Truncated") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
ggsave("mjm_scenII_trunc_coll_condi.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)


ratio_dat <- data.frame(
  it = rep(rep(start:stop, times = nmarker), 2),
  vals = c(ratio_x, ratio_i),
  type = factor(rep(c(0, 1), each = nmarker * length(start:stop)),
                labels = c("XTX", "I_S")),
  ev = factor(seq_len(nmarker), labels = paste("EV", seq_len(nmarker)))
)
ggplot(ratio_dat, aes(x = it, y = vals, col = type,
                      group = interaction(ev, type))) +
  geom_line() +
  theme_bw() +
  ggtitle("Ratio of Eigenvalues for SimScen II It 100 Truncated") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)
ggsave("mjm_scenII_trunc_coll_ratio.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)


# Calcualate Model With True FPCs (E) --------------------------------------

# Load the data
load(file.path(server_wd, "JMbamlss/simulation", dat_setting, 
               "data/d100.Rdata"))

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
m1sp1 <- splinefun(x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 
                        0.5))
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

# Use the true data generating FPCs
# Use all FPCs 
nfpc <- nObs(d_sim$fpc_base)
mfpca_est_list <- lapply(seq_len(nfpc), function (i) {
  list(functions = extractObs(m$functions, i),
       values = m$values[i])
})

# Prepare objects for model fit
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker +",
    " obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], scale = 'FALSE'))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, 
                data = d_sim$data, 
                timevar = "obstime", maxit = 1500, sampler = FALSE,
                coll = TRUE, verbose = TRUE)


# Analyze Collinearity Model True/No Trunc --------------------------------


nmarker <- length(m$functions)
start <- max(which(sapply(b_est$model.stats$optimizer$coll, is.null))) + 1
stop <- length(b_est$model.stats$optimizer$coll)

# Use data matrix
col_x <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  X <- matrix(x$X, ncol = nmarker)
  svd(crossprod(X))$d
})
condi_x <- apply(col_x, 2, function (x) x[1] / x[nmarker])
ratio_x <- apply(col_x, 2, function (x) x / sum(x))

# Use Hesse matrix
col_i <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  H <- matrix(x$I, ncol = nmarker)
  s <- diag(diag(H)^(-1/2))
  svd(s %*% H %*% s)$d
})
condi_i <- apply(col_i, 2, function (x) x[1] / x[nmarker])
ratio_i <- apply(col_i, 2, function (x) x / sum(x))

condi_dat <- data.frame(
  it = rep(start:stop, times = 2),
  vals = c(condi_x, condi_i),
  type = factor(rep(c(0, 1), each = length(start:stop)), labels = c("XTX", "I_S"))
)
ggplot(condi_dat, aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Condition Numbers for SimScen II It 100 True") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
ggsave("mjm_scenII_true_coll_condi.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)


ratio_dat <- data.frame(
  it = rep(rep(start:stop, times = nmarker), 2),
  vals = c(ratio_x, ratio_i),
  type = factor(rep(c(0, 1), each = nmarker * length(start:stop)),
                labels = c("XTX", "I_S")),
  ev = factor(seq_len(nmarker), labels = paste("EV", seq_len(nmarker)))
)
ggplot(ratio_dat, aes(x = it, y = vals, col = type,
                      group = interaction(ev, type))) +
  geom_line() +
  theme_bw() +
  ggtitle("Ratio of Eigenvalues for SimScen II It 100 True") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)
ggsave("mjm_scenII_true_coll_ratio.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)
