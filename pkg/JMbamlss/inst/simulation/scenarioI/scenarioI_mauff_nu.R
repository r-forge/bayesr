
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
  "workstation" = paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.",
                         "hu-berlin.de,share=volkmana.hub/JMbamlss/simulation/"),
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
  
  d_rirs <- pivot_longer(dat, y1:y6, names_to = "marker", values_to = "y") %>%
    mutate(marker = factor(marker, labels = paste0("m", 1:6)),
           id = factor(id)) %>%
    arrange(marker, id, year) %>%
    as.data.frame()
  
  try_obj <- try({
    
    # True FPCs ---------------------------------------------------------------
    
    D_num <- matrix(NA, nrow = 12, ncol = 12)
    D_num[1:4, 1:4] <- matrix(c(0.919,  0.003,  0.220, -0.008,
                                0.003,  0.006,  0.005, -0.002,
                                0.220,  0.005,  0.421, -0.013,
                                -0.008, -0.002, -0.013,  0.008), byrow = TRUE)
    D_num[5:8, 5:8] <- matrix(c(0.551,  0.007, -0.141,  0.015,
                                0.007,  0.014, -0.005, -0.005,
                                -0.141, -0.005,  0.204, -0.042,
                                0.015, -0.005, -0.042,  0.043), byrow = TRUE)  
    D_num[9:12, 9:12] <- matrix(c(0.110, -0.028,  0.176,  0.021,
                                  -0.028,  0.035, -0.026, -0.003,
                                  0.176, -0.026,  3.580,  0.040,
                                  0.021, -0.003,  0.040,  0.197), byrow = TRUE)  
    D_num[1:4, 5:8] <- matrix(c(0.308,  0.005,  0.265, -0.013,
                                0.012,  0.001,  0.012,  0.002,
                                -0.073, -0.007, -0.172,  0.018,
                                0.007,  0.005,  0.012, -0.013), byrow = TRUE)
    D_num[1:4, 9:12] <- matrix(c(0.049,  0.005,  0.124, -0.013,
                                 -0.007, -0.005, -0.010,  0.012,
                                 0.698,  0.006,  0.680, -0.027,
                                 0.056,  0.006,  0.034,  0.004), byrow = TRUE)  
    D_num[5:8, 9:12] <- matrix(c(0.095,  0.004, -0.144,  0.032,
                                 -0.013,  0.005,  0.035, -0.037,
                                 0.826,  0.018, -0.263,  0.024,
                                 0.077,  0.019, -0.015,  0.006), byrow = TRUE)
    D_num[5:8, 1:4] <- t(D_num[1:4, 5:8])
    D_num[9:12, 1:4] <- t(D_num[1:4, 9:12])
    D_num[9:12, 5:8] <- t(D_num[5:8, 9:12])
    
    D_numsym <- forceSymmetric(D_num) 
    D <- D_numsym
    
    # Basis functions on each dimension
    seq1 <- seq(0, max(d_rirs$Time), length.out = 101)
    b_funs <- rep(list(funData(argvals = seq1,
                               X = matrix(c(rep(1, length(seq1)), seq1),
                                          byrow = TRUE, ncol = length(seq1)))), 6)
    
    mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs)
    
    
    
    # TRUE MFPC Model ---------------------------------------------------------
    
    mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                         function (i, mfpca = mfpca_tru) {
                           list(functions = extractObs(mfpca$functions, i),
                                values = mfpca$values[i])
                         })
    
    mfpca_list_w <- lapply(seq_along(mfpca_tru$values), 
                           function (i, mfpca = mfpca_tru) {
                             list(functions = extractObs(mfpca$functions, i) *
                                    sqrt(mfpca$values[[i]]),
                                  values = mfpca$values[i])
                           })
    
    # Prepare objects for model fit
    d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs, n = 12,
                                         obstime = "year")
    d_rirs_tru_w <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs, n = 12,
                                         obstime = "year", eval_weight = TRUE)
    
    nfpc <- 12
    f_tru <- list(
      Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
      gamma ~ 1 + group,
      as.formula(paste0(
        "mu ~ -1 + marker + year:marker +",
        paste0(lapply(seq_along(1:5), function(x) {
          paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
                 "mfpca_list[[", x, "]]))")
        }), collapse = " + "))),
      sigma ~ -1 + marker,
      alpha ~ -1 + marker
    )
    f_prio <- list(
      Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
      gamma ~ 1 + group,
      as.formula(paste0(
        "mu ~ -1 + marker + year:marker +",
        paste0(lapply(seq_along(1:5), function(x) {
          paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
                 "mfpca_list[[", x, "]], a = 1.0001, b = 1/mfpca_list[[", x,
                 "]]$values))")
        }), collapse = " + "))),
      sigma ~ -1 + marker,
      alpha ~ -1 + marker
    )
    f_tru_w <- list(
      Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
      gamma ~ 1 + group,
      as.formula(paste0(
        "mu ~ -1 + marker + year:marker +",
        paste0(lapply(seq_along(1:5), function(x) {
          paste0("s(id, wfpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
                 "mfpca_list[[", x, "]]))")
        }), collapse = " + "))),
      sigma ~ -1 + marker,
      alpha ~ -1 + marker
    )
    
    # Model fit
    t_est <- system.time(
      b_prio <- bamlss(f_prio, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                      timevar = "year", maxit = 1500, n.iter = 1500, nu = 1,
                      burnin = 500, thin = 5, update_nu = TRUE, verbose  = TRUE)
    )
    attr(b_est, "comp_time") <- t_est
    attr(b_est, "FPCs") <- mfpca_tru
    attr(b_est, "nfpc") <- nfpc
    save(b_est, file = paste0(server_wd, setting, "/bamlss_tru_nu/b", i, 
                              ".Rdata"))
  }, silent = TRUE)
  if(class(try_obj) == "try-error") try_obj else i
  
  b_est_w <- bamlss(f_tru_w, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru_w, 
                  timevar = "year", maxit = 1500, n.iter = 5500,
                  burnin = 500, thin = 5, update_nu = TRUE, verbose = TRUE)
}



# Simulation --------------------------------------------------------------

cat(paste("Parallelization:", start, "-", stop, "on", number_cores, "cores.\n"))
simulation <- mclapply(start:stop, parallel_bamlss_est, mc.cores = number_cores)
