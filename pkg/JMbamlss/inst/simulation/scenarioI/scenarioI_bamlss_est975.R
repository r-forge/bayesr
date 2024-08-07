
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
start <- 100
stop <- 299
number_cores <- 2
setting <- "scen_I_230719"
Sys.time()
sessionInfo()


# Simulation function -----------------------------------------------------

parallel_bamlss_est <- function(i) {
  set.seed(i)
  
  # Load the data
  d_rirs <- readRDS(file.path(results_wd, setting, "data",
                              paste0("d", i, ".rds")))

  try_obj <- try({
    
    # Estimate the model using estimated FPCs
    few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                     function (x) any(x < 2))
    long_obs <- d_rirs$data %>%
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
      d_rirs$data %>%
        filter(id %in% take) %>% 
        droplevels(), 
      uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
      npc = 2, nbasis = 4
    )
    vals <- which(mfpca_est$values > 0)
    nfpc <- min(which(
      cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
    mfpca_est_list <- lapply(vals[seq_len(nfpc)], 
                             function (i, mfpca = mfpca_est) {
      list(functions = extractObs(mfpca$functions, i),
           values = mfpca$values[i])
    })
    
    # Prepare objects for model fit
    d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
    f_est <- list(
      Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
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
    t_est <- system.time(
      b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                      timevar = "obstime", maxit = 1500, n.iter = 5500,
                      burnin = 500, thin = 5)
    )
    attr(b_est, "comp_time") <- t_est
    attr(b_est, "FPCs") <- mfpca_est
    attr(b_est, "nfpc") <- nfpc
    saveRDS(b_est, file = file.path(results_wd, setting, "bamlss_est975", 
                                    paste0("b", i, ".rds")))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
