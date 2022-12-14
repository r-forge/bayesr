
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
server_wd <- switch(location,
                    "workstation" = paste0("/run/user/1000/gvfs/smb-share:",
                                           "server=clapton.wiwi.hu-berlin.de,",
                                           "share=volkmana.hub/JMbamlss/",
                                           "simulation/"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation/",
                    "server_windows" = "H:/JMbamlss/simulation/")


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
start <- 1
stop <- 49
number_cores <- 2
setting <- "scen_mauff"
Sys.time()
sessionInfo()


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))
  
  # Get the quantile-based knots for comparability
  kn <- mgcv::smoothCon(mgcv::s(Time, k = 20, bs = "ps"), 
                        data = dat)[[1]]$knots
  
  # Estimate the model using JMbayes
  try_obj <- try({
    CoxFit <- coxph(Surv(Time, event) ~ group, 
                    data = dat.id)
    lm1 <- lme(y1 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm2 <- lme(y2 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm3 <- lme(y3 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm4 <- lme(y4 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm5 <- lme(y5 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    lm6 <- lme(y6 ~ year, random = ~ year | id, 
               data = dat, na.action = na.omit,
               control = lmeControl(opt = "optim"))
    
    t_jm <- system.time(
      jmb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "year",
                n_iter = 1900L, n_burnin = 1000L, n_thin = 3L, 
                cores = 1, n_chains = 1, 
                GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                save_random_effects = TRUE)
    )
    attr(jmb, "comp_time") <- t_jm
    save(jmb, file = paste0(server_wd, setting, "/jmb/jmb", i, 
                              ".Rdata"))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
