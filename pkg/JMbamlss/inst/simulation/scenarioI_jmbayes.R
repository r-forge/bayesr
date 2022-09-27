
# Set Up Small Simulation -------------------------------------------------



# Set up R session --------------------------------------------------------



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
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")


# Setting for the simulation
start <- 100
stop <- 149
number_cores <- 10
setting <- "scen_I_130922"
Sys.time()
sessionInfo()


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  load(paste0("simulation/", setting, "/data/d", i, ".Rdata"))
  
  # Estimate the model using JMbayes
  d_rirs_jmb <- d_rirs$data %>% 
    pivot_wider(names_from = marker, values_from = y)
  CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                  data = d_rirs$data_short %>% filter(marker == "m1"))
  lm1 <- lme(m1 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm2 <- lme(m2 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm3 <- lme(m3 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm4 <- lme(m4 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm5 <- lme(m5 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  lm6 <- lme(m6 ~ obstime * x3, random = ~ obstime | id, 
             data = d_rirs_jmb, na.action = na.omit,
             control = lmeControl(opt = "optim"))
  
  t_jm <- system.time(
    jmb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "obstime",
              n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
              control = list(cores = 1, n_chains = 1))
  )
  attr(jmb, "comp_time") <- t_jm
  save(jmb, file = paste0("simulation/", setting, "/jmb/jmb_", i, 
                          ".Rdata"))
  NULL
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
