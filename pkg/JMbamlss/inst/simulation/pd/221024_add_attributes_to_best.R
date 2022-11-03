
# Set Up Small Simulation -------------------------------------------------



# Set up R session --------------------------------------------------------



# Specify location
location <- "server_linux"
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
m_1309 <- list.files("simulation/scen_I_130922/bamlss_est/")
m1_2210 <- list.files("simulation/scen_I_221019/bamlss_est_1/")
m_2210 <- list.files("simulation/scen_I_221019/bamlss_est/")


# Simulation function -----------------------------------------------------

add_attr_bamlss_est <- function(i, setting, folder, d_set = "scen_I_130922") {

  # Load the data
  load(paste0("simulation/", d_set, "/data/d", substr(i, 2, 10)))
  load(paste0("simulation/", setting, "/", folder, "/", i))
  

    # Estimate the model using estimated FPCs
    few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                     function (x) any(x < 2))
    mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                                 filter(id %in% paste(which(!few_obs))) %>% 
                                 droplevels(), 
                               uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                               npc = 2, nbasis = 4)
    attr(b_est, "FPCs") <- mfpca_est
    save(b_est, file = paste0("simulation/", setting, "/", folder, "/", i))
    
    cat(i, "\n")
}

vadd_attr_bamlss_est <- Vectorize(add_attr_bamlss_est, vectorize.args = "i")


vadd_attr_bamlss_est(m_1309, setting = "scen_I_130922", folder = "bamlss_est")
vadd_attr_bamlss_est(m1_2210, setting = "scen_I_221019", 
                     folder = "bamlss_est_1")
vadd_attr_bamlss_est(m_2210, setting = "scen_I_221019", folder = "bamlss_est")
