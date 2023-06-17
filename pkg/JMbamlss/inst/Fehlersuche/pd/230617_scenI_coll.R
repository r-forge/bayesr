
# Collinearity in SimScen I -----------------------------------------------


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
dat_setting <- "scen_I_130922"


# Always
library(survival)
library(JMbayes2)
library(tidyverse)
devtools::load_all(".")



# Calculate Iteration 100 -------------------------------------------------

# Load the data
load(file.path(server_wd, "JMbamlss/simulation", dat_setting, 
               "data/d100.Rdata"))

# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                 function (x) any(x < 2))
mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                             filter(id %in% paste(which(!few_obs))) %>% 
                             droplevels(), 
                           uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                           npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
d_rirs_est <- attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
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
b_est <- bamlss(f_est, family = mjm_bamlss, data = d_rirs_est, 
                timevar = "obstime", maxit = 1500, sampler = FALSE,
                coll = TRUE, verbose = TRUE)



# Analyze Collinearity Model 100 ------------------------------------------


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
  ggtitle("Condition Numbers for SimScen I Iteration 100") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)


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
  ggtitle("Ratio of Eigenvalues for SimScen I Iteration 100") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)


# Calculate Iteration 151 -------------------------------------------------

# Load the data
load(file.path(server_wd, "JMbamlss/simulation", dat_setting, 
               "data/d151.Rdata"))

# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                 function (x) any(x < 2))
mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                             filter(id %in% paste(which(!few_obs))) %>% 
                             droplevels(), 
                           uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                           npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
d_rirs_est <- attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
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
b_est <- bamlss(f_est, family = mjm_bamlss, data = d_rirs_est, 
                timevar = "obstime", maxit = 1500, sampler = FALSE,
                coll = TRUE, verbose = TRUE)



# Analyze Collinearity Model 151 ------------------------------------------


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
  ggtitle("Condition Numbers for SimScen I Iteration 151") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)


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
  ggtitle("Ratio of Eigenvalues for SimScen I Iteration 151") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)

