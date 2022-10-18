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

# ggplot(pbc2, aes(x = year, y = log(serBilir), group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = log(serChol), group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = albumin, group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = log(alkaline), group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = log(SGOT), group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = log(platelets), group = id)) +
#   geom_line()
# ggplot(pbc2, aes(x = year, y = log(prothrombin), group = id)) +
#   geom_line()

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
                           uni_mean = "y ~ 1 + s(obstime) + sex")#,
#npc = 2, nbasis = 4)
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
  as.formula(paste0(
    "y ~ -1 + marker + s(obstime, by = marker) + sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker
)

# Model fit
b_est <- bamlss(f_est, data = p_long, 
                maxit = 1500, n.iter = 5500,
                burnin = 500, thin = 5)
