

# Simulation Scenario II - 22/12/09 ---------------------------------------




# Set up R session --------------------------------------------------------


# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if(location == "server_linux") "./simulation"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation")
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
setting <- "scen_II_230719"
Sys.time()
sessionInfo()
i <- 100


# Simulation function -----------------------------------------------------


set.seed(i)

# Load the data
d_sim <- readRDS(file.path(results_wd, setting, "data",
                           paste0("d", i, ".rds")))


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
  d_sim$data %>%
    filter(id %in% take) %>% 
    droplevels(), 
  uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
  npc = 3, nbasis = 10
)
vals <- which(mfpca_est$values > 0)

# Use all FPCs 
nfpc_est1 <- length(mfpca_est$values)
nfpc_est975 <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
nfpc_est95 <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
mfpca_est_list <- lapply(vals[seq_len(nfpc_est1)], function (i,
                                                        mfpca = mfpca_est){
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Remove old FPC variables
d_sim_est <- d_sim$data[, -grep("fpc\\.", names(d_sim$data))]



# Full model --------------------------------------------------------------


# Prepare objects for model fit
d_sim_est1 <- JMbamlss:::attach_wfpc(mfpca_est, d_sim_est, n = nfpc_est1)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc_est1), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est1 <- system.time(
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_sim_est1,
                  timevar = "obstime", maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5)
)
rm(b_est)



# Truncated Model (975) ---------------------------------------------------

# Prepare objects for model fit
d_sim_est975 <- JMbamlss:::attach_wfpc(mfpca_est, d_sim_est, n = nfpc_est975)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc_est975), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est975 <- system.time(
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_sim_est975,
                  timevar = "obstime", maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5)
)
rm(b_est)


# Truncated Model (95) ---------------------------------------------------

# Prepare objects for model fit
d_sim_est95 <- JMbamlss:::attach_wfpc(mfpca_est, d_sim_est, n = nfpc_est95)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc_est95), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est95 <- system.time(
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_sim_est95,
                  timevar = "obstime", maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5)
)
rm(b_est)



# JMbayes2 model ----------------------------------------------------------

# Get the quantile-based knots for comparability
kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                      data = d_sim$data)[[1]]$knots

# Estimate the model using JMbayes
d_sim_jmb <- d_sim$data %>% 
  pivot_wider(names_from = marker, values_from = y)
t_jm <- system.time({
  CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                  data = d_sim$data_short %>% filter(marker == "m1"))
  lm1 <- lme(m1 ~ obstime * x3, 
             random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
             data = d_sim_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm2 <- lme(m2 ~ obstime * x3, 
             random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
             data = d_sim_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  jmb <- jm(CoxFit, list(lm1, lm2), time_var = "obstime",
            n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
            cores = 1, n_chains = 1, 
            GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
            save_random_effects = TRUE)
})


# > t_est1
# user   system  elapsed 
# 2581.787   36.698 2617.837 
# > t_est975
# user   system  elapsed 
# 1747.960    6.238 1753.776 
# > t_est95
# user   system  elapsed 
# 1374.865    0.503 1375.025 
# > t_jm
# user  system elapsed 
# 57.993   0.144  58.142 