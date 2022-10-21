
# PBC Data Analysis -------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
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
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/MJM_opt.R")
source("R/opt_updating_cpp.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing_cpp.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/fpca.R")
source("R/mfpca_sim.R")
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")

# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years, obstime = year, marker = factor(marker)) %>%
  select(id, survtime, event, sex, age, marker, obstime, y) %>%
  na.omit() %>%
  as.data.frame()

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis)

# Estimate the model using estimated FPCs
few_obs <- apply(table(p_long$id, p_long$marker), 1,
                 function (x) any(x < 2))
mfpca_est <- preproc_MFPCA(p_long %>%
                             filter(id %in% paste(which(!few_obs))) %>% 
                             droplevels(), 
                           uni_mean = "y ~ 1 + s(obstime) + sex")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long <- attach_wfpc(mfpca_est, p_long, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker) + sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est <- system.time(
  b_est <- bamlss(f_est, family = mjm_bamlss, data = p_long,
                  timevar = "obstime", maxit = 15,
                  sampler = FALSE, verbose = TRUE)
)
# No difference in time if verbose = FALSE
# t_est_noverb <- system.time(
#   b_est <- bamlss(f_est, family = mjm_bamlss, data = p_long,
#                   timevar = "obstime", maxit = 15,
#                   sampler = FALSE, verbose = FALSE)
# )
# 

source("R/opt_updating_cpp_new.R")
t_new <- system.time(
  b_new <- bamlss(f_est, family = mjm_bamlss, data = p_long,
                  timevar = "obstime", maxit = 15,
                  sampler = FALSE, verbose = TRUE)
)
