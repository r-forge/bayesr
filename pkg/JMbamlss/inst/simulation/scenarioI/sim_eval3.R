

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


debugonce(JMbamlss:::sim_jmbamlss_eval)
m_tru <- JMbamlss:::sim_jmbamlss_eval(
  wd = file.path(server_wd, "scen_I_230719/"), model_wd = "bamlss_tru/",
  data_wd = "data/", name = "bamlss_tru")

debugonce(JMbamlss:::sim_jmbayes_eval)
m_jmb <- JMbamlss:::sim_jmbayes_eval(
  wd = file.path(server_wd, "scen_I_230719/"), model_wd = "jmb/",
  data_wd = "data/", name = "jmb")

m_tru <- list.files(path = paste0(wd, "/bamlss_tru"))
m_est1 <- list.files(path = paste0(wd, "/bamlss_est1"))
m_est95 <- list.files(path = paste0(wd, "/bamlss_est95"))


# load(paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
#             "simulation/scen_I_130922/bamlss_tru/b101.Rdata"))
# mfpca_est <- attr(b_est, "FPCs")

res_tru <- lapply(m_tru, function (x) {
  
  # Load the data set and extract information about it
  m <- readRDS(paste0(wd, "/bamlss_tru/", x))
  load(paste0(wd, "/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  #mfpca_est <- attr(b_est, "FPCs")
  vals <- which(mfpca_est$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
  newdat <- attach_wfpc(mfpca_est, d_rirs$data_full, n = nfpc)
  
  list("lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.025),
                             "Mean" = rowMeans(mcmc_lambga),
                             "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                             probs = 0.975)),
       "alpha" = cbind(predict(b_est, model = "alpha", FUN = c95)[marks, ],
                       data.frame("marker" = b_est$model.frame$marker[marks])),
       "mu" = cbind(predict(b_est, model = "mu", FUN = c95), 
                    data.frame("marker" = b_est$model.frame$marker)),
       "sigma" = cbind(predict(b_est, model = "sigma", FUN = c95)[marks, ],
                       data.frame("marker" = b_est$model.frame$marker[marks])),
       "mu_long" = cbind(predict(b_est, model = "mu", FUN = c95, 
                                 newdata = newdat),
                         data.frame("marker" = newdat$marker,
                                    "obstime" = newdat$obstime)))
})
save(res_full, file = paste0(wd, "/res_full.Rdata"))


res_975 <- lapply(m_975, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(wd, "/bamlss_tru/", x))
  load(paste0(wd, "/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  #mfpca_est <- attr(b_est, "FPCs")
  vals <- which(mfpca_est$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
  newdat <- attach_wfpc(mfpca_est, d_rirs$data_full, n = nfpc)
  
  list("lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.025),
                             "Mean" = rowMeans(mcmc_lambga),
                             "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                             probs = 0.975)),
       "alpha" = cbind(predict(b_est, model = "alpha", FUN = c95)[marks, ],
                       data.frame("marker" = b_est$model.frame$marker[marks])),
       "mu" = cbind(predict(b_est, model = "mu", FUN = c95), 
                    data.frame("marker" = b_est$model.frame$marker)),
       "sigma" = cbind(predict(b_est, model = "sigma", FUN = c95)[marks, ],
                       data.frame("marker" = b_est$model.frame$marker[marks])),
       "mu_long" = cbind(predict(b_est, model = "mu", FUN = c95, 
                                 newdata = newdat),
                         data.frame("marker" = newdat$marker,
                                    "obstime" = newdat$obstime)))
})
save(res_975, file = paste0(wd, "/res_975.Rdata"))



d_full <- paste0("d", substr(m_full, 2, 10))
d_975 <- paste0("d", substr(m_975, 2, 10))

sim_dat_full <- lapply(d_975, function (x) {
  
  load(paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
                      "simulation/scen_I_130922/", "data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sim_dat_full, file = paste0(wd, "/sim_dat_full.Rdata"))

sim_dat_975 <- lapply(d_full, function (x) {
  
  load(paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
              "simulation/scen_I_130922/", "data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sim_dat_975, file = paste0(wd, "/sim_dat_975.Rdata"))

load(paste0(wd, "/sim_dat_full.Rdata"))
load(paste0(wd, "/sim_dat_975.Rdata"))
load(paste0(wd, "/res_full.Rdata"))
load(paste0(wd, "/res_975.Rdata"))

r_full <- sim_results(res_full, sim_dat_full, name = "full")
r_975 <- sim_results(res_975, sim_dat_975, name = "975")
save(r_full, r_975, file = paste0(wd, "/results.Rdata"))
