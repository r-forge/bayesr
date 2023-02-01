

# Simulation Scenario II - 22/12/09 ---------------------------------------




# Set up R session --------------------------------------------------------


# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/fpc_estimation"
        else "H:/fpc_estimation")
  results_wd <- if(location == "server_linux") "./simulation/"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/fpc_estimation/simulatio",
                       "n/")
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
number_cores <- 5
setting <- "230124_parallel_to_simscenII/"
dat_setting <- "230124_parallel_to_simscenII/"
model_type <- "A_uncens/"
dir.create(file.path(results_wd, setting, model_type), showWarnings = FALSE)
Sys.time()
sessionInfo()


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_sim <- readRDS(paste0(results_wd, dat_setting, "/data/d", i, ".rds"))
  
  try_obj <- try({
    
    # Estimate the model using estimated FPCs
    # remove observation with less than 3 longitudinal observations
    few_obs <- apply(table(d_sim$data_uncens$id, d_sim$data_uncens$marker), 1, 
                     function (x) any(x < 3))
    long_obs <- d_sim$data_uncens %>%
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
      d_sim$data_uncens %>% filter(id %in% take) %>% droplevels(), 
      uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
      npc = 3, nbasis = 10, fit = TRUE, save_uniGAM = TRUE)
    
    saveRDS(mfpca_est, file = paste0(results_wd, setting, model_type, "m", i, 
                                     ".rds"))
    
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, 
          "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est,
                       mc.cores = number_cores)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# Now estimation of the randomly censored data ----------------------------

model_type <- "A_cens/"
dir.create(file.path(results_wd, setting, model_type), showWarnings = FALSE)


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_sim <- readRDS(paste0(results_wd, dat_setting, "/data/d", i, ".rds"))
  
  try_obj <- try({
    
    # Estimate the model using estimated FPCs
    # remove observation with less than 3 longitudinal observations
    few_obs <- apply(table(d_sim$data_cens$id, d_sim$data_cens$marker), 1, 
                     function (x) any(x < 3))
    long_obs <- d_sim$data_cens %>%
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
      d_sim$data_cens %>% filter(id %in% take) %>% droplevels(), 
      uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
      npc = 3, nbasis = 10, fit = TRUE, save_uniGAM = TRUE)
    
    saveRDS(mfpca_est, file = paste0(results_wd, setting, model_type, "m", i, 
                                     ".rds"))
    
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, 
          "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est,
                       mc.cores = number_cores)



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# Now estimation of the informative censoring -----------------------------

model_type <- "A_jm/"
dat_setting <- "scen_II_230117/"
dir.create(file.path(results_wd, setting, model_type), showWarnings = FALSE)


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_sim <- readRDS(paste0(results_wd, "../../JMbamlss/simulation/",
                          dat_setting, "/data/d", i, ".rds"))
  
  try_obj <- try({
    
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
      npc = 3, nbasis = 10, fit = TRUE, save_uniGAM = TRUE)
    
    saveRDS(mfpca_est, file = paste0(results_wd, setting, model_type, "m", i, 
                                     ".rds"))
    
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, 
          "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est,
                       mc.cores = number_cores)
