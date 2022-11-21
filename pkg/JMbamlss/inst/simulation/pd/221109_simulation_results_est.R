location <- "workstation"

if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if (location == "server_linux") {
    paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/simulation/")
  } else NULL
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
}

library(survival)
library(MFPCA)
library(tidyverse)
library(devtools)
load_all()



# Prediction Part ---------------------------------------------------------

m_est_975 <- list.files(path = paste0(results_wd, "/scen_I_221019/bamlss_est"))
m_est_1 <- list.files(path = paste0(results_wd, "/scen_I_221019/bamlss_est_1"))
m_jmb <- list.files(path = paste0(results_wd, "/scen_I_221019/jmb"))

pred_est_975 <- sim_bamlss_predict(m_est_975, results_wd, 
                                   m_setting = "/scen_I_221019/", 
                                   d_setting = "/scen_I_130922/", 
                                   folder = "bamlss_est/")
pred_est_1 <- sim_bamlss_predict(m_est_1, results_wd, 
                                 m_setting = "/scen_I_221019/", 
                                 d_setting = "/scen_I_130922/", 
                                 folder = "bamlss_est_1/")
pred_jmb <- sim_jmb_predict(m_jmb, results_wd,
                            m_setting = "/scen_I_221019/",
                            d_setting = "/scen_I_130922/",
                            folder = "jmb/")



# Simulated Data Part -----------------------------------------------------

d_est_975 <- paste0("d", substr(m_est_975, 2, 10))
d_est_1 <- paste0("d", substr(m_est_1, 2, 10))
d_jmb <- paste0("d", substr(m_jmb, 5, 13))

sdat_est_975 <- lapply(d_est_975, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  # Remove mu_long times that are higher than maximum observation time
  max_t <- sapply(split(d_rirs$data, d_rirs$data$marker),
                  function(x) max(as.numeric(names(table(x$obstime)))))
  if (any(max_t < 1)) {
    for(m in seq_along(max_t)) {
      d_rirs$data_full <- subset(d_rirs$data_full, 
                                 !(d_rirs$data_full$marker == names(max_t)[m] &
                                     d_rirs$data_full$obstime > max_t[m]))
    }
  }
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_est_975, file = paste0(results_wd, "/sdat_est_975.Rdata"))

sdat_est_1 <- lapply(d_est_1, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  # Remove mu_long times that are higher than maximum observation time
  max_t <- sapply(split(d_rirs$data, d_rirs$data$marker),
                  function(x) max(as.numeric(names(table(x$obstime)))))
  if (any(max_t < 1)) {
    for(m in seq_along(max_t)) {
      d_rirs$data_full <- subset(d_rirs$data_full, 
                                 !(d_rirs$data_full$marker == names(max_t)[m] &
                                     d_rirs$data_full$obstime > max_t[m]))
    }
  }
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_est_1, file = paste0(results_wd, "/sdat_est_1.Rdata"))


sdat_jmb <- lapply(d_jmb, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_jmb, file = paste0(results_wd, "/sdat_jmb.Rdata"))



# Comparing Predictions to Data -------------------------------------------

r_est_975 <- sim_results(pred_est_975, sdat_est_975, name = "est_975")
r_est_1 <- sim_results(pred_est_1, sdat_est_1, name = "est_1")
r_jmb <- sim_results(pred_jmb, sdat_jmb, name = "jmb")
save(r_est_975, r_est_1, r_jmb, 
     file = paste0(results_wd, "/results_est.Rdata"))

