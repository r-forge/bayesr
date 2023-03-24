
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
  d_rirs <- readRDS(paste0(server_wd, setting, "/data/data", i, ".rds"))
  
  
  try_obj <- try({
    
    # Estimate the model using estimated FPCs
    few_obs <- apply(table(d_rirs$id, d_rirs$marker), 1, 
                     function (x) any(x < 2))
    long_obs <- d_rirs %>%
      group_by(id, marker) %>%
      summarize(maxobs = max(year), .groups = "drop_last") %>%
      ungroup(marker) %>%
      summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
      ungroup() %>%
      filter(minmaxobs > 0.1) %>%
      select(id) %>%
      unlist() %>%
      paste()
    take <- intersect(long_obs, paste(which(!few_obs)))
    
    mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs %>%
                                 filter(id %in% take) %>% 
                                 droplevels(), 
                                 time = "year", weights = TRUE,
                               uni_mean = "y ~ 1 + year",
                               npc = 2, nbasis = 4)
    vals <- which(mfpca_est$values > 0)
    
    uni_norms <- lapply(mfpca_est$functions, norm)
    nfpc <- max(sapply(uni_norms, function (n) {
      min(which(
        cumsum(mfpca_est$values[vals]*n)/sum(mfpca_est$values[vals]*n) > 0.9
      ))
    }))
    mfpca_est_list <- lapply(vals[seq_len(nfpc)], 
                             function (i, mfpca = mfpca_est) {
      list(functions = extractObs(mfpca$functions, i),
           values = mfpca$values[i])
    })
    
    # Prepare objects for model fit
    d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs, n = nfpc,
                                         obstime = "year")
    f_est <- list(
      Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
      gamma ~ 1 + group,
      as.formula(paste0(
        "mu ~ -1 + marker + year:marker +",
        paste0(lapply(seq_len(nfpc), function(x) {
          paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
                 "mfpca_est_list[[", x, "]]))")
        }), collapse = " + "))),
      sigma ~ -1 + marker,
      alpha ~ -1 + marker
    )
    
    # Model fit
    t_est <- system.time(
      b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                      timevar = "year", maxit = 1500, n.iter = 900,
                      burnin = 1000, thin = 3, verbose = TRUE, nu = 1)
    )
    attr(b_est, "comp_time") <- t_est
    attr(b_est, "FPCs") <- mfpca_est
    attr(b_est, "nfpc") <- nfpc
    save(b_est, file = paste0(server_wd, setting, "/bamlss_est/b", i, 
                              ".Rdata"))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
