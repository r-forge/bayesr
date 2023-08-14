
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
                                           "simulation"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation",
                    "server_windows" = "H:/JMbamlss/simulation")


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



# Extract FPCs from models ------------------------------------------------


extract_fpc_info_i <- function (i, results_wd, setting, folder) {
  
  b_est <- readRDS(file.path(results_wd, setting, folder, i))
  
  list("nfpc" = attr(b_est, "nfpc"),
       "mfpca" = attr(b_est, "FPCs"))
}
extract_fpc_info <- Vectorize(extract_fpc_info_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
fpc_info <- function (results_wd, setting, folder) {
  
  m <- list.files(file.path(results_wd, setting, folder))
  extract_fpc_info(m, results_wd = results_wd, setting = setting,
                   folder = folder)
  
}

debug(fpc_info)

simI <- fpc_info(results_wd = server_wd,
                 setting = "scen_I_230719",
                 folder = "bamlss_est95")
simII <- fpc_info(results_wd = server_wd,
                  setting = "scen_II_230719",
                  folder = "bamlss_est95")



# Number of FPCs in the Model ---------------------------------------------

summary(sapply(simI, "[[", "nfpc"))
summary(sapply(simII, "[[", "nfpc"))


# Extract Eigenfunctions and Combine to Data.frame ------------------------

simI_tru <- extract_fpc_info_i("b100.rds", server_wd, "scen_I_230719", 
                               "bamlss_tru")
simII_tru <- extract_fpc_info_i("b100.rds", server_wd, "scen_II_230719", 
                                "bamlss_tru")


eigenfun_norm_i <- function(mfpca, mfpca_tru) {
  sapply(seq_along(mfpca$functions), function (fu_it) {
    tr <- extractObs(mfpca_tru$mfpca$functions, fu_it)
    es <- flipFuns(tr, extractObs(mfpca$functions, fu_it))
    norm(tr - es)
  })
}
  dat <- lapply(mfpca$functions@.Data, function(m) {
    data.frame(t = m@argvals[[1]], 
               X = c(t(m@X)), 
               fpc = factor(paste0("fpc.", rep(seq_along(mfpca$values), 
                                               each = length(m@argvals[[1]]))),
                            levels = paste0("fpc.", seq_along(mfpca$values))))
  })
  do.call(rbind, Map(cbind, marker = factor(paste0("m", seq_along(dat))), dat))
}

eigenfun_data_list <- Vectorize(eigenfun_data_i, SIMPLIFY = FALSE)

eigenfun_data <- function(mfpca_list) {
  dat <- eigenfun_data_list(mfpca_list)
  do.call(rbind, Map(cbind, it = factor(seq_along(dat)), dat))
}