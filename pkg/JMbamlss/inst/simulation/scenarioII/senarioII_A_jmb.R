

# Simulation Scenario II - 22/12/09 ---------------------------------------


# JMbayes estimation


# Set up R session --------------------------------------------------------


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



# Setting for the simulation
start <- 100
stop <- 199
number_cores <- 2
setting <- "scen_II_221209/"
dat_setting <- "scen_II_221209/"
model_type <- "A_jmb/"
Sys.time()
sessionInfo()



# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  load(paste0(results_wd, dat_setting, "data/d", i, ".Rdata"))
  
  # Get the quantile-based knots for comparability
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                        data = d_rirs$data)[[1]]$knots
  
  # Estimate the model using JMbayes
  d_rirs_jmb <- d_rirs$data %>% 
    pivot_wider(names_from = marker, values_from = y)
  try_obj <- try({
    CoxFit <- coxph(Surv(survtime, event) ~ x3, 
                    data = d_rirs$data_short %>% filter(marker == "m1"))
    lm1 <- lme(m1 ~ obstime * x3, 
               random = ~ bs(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
               data = d_rirs_jmb, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm2 <- lme(m2 ~ obstime * x3, 
               random = ~ ns(obstime, df = 3, Boundary.knots = c(0, 1)) | id, 
               data = d_rirs_jmb, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    
    t_jm <- system.time(
      jmb <- jm(CoxFit, list(lm1, lm2), time_var = "obstime",
                n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
                cores = 1, n_chains = 1, 
                GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                save_random_effects = TRUE)
    )
    attr(jmb, "comp_time") <- t_jm
    save(jmb, file = paste0(results_wd, setting, model_type, "jmb", i, 
                            ".Rdata"))
      
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
