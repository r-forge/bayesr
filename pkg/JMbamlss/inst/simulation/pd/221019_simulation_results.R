location <- "workstation"


if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


if (location == "server_linux") {
  results_wd <- paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
                       "simulation/")
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
source("R/sim_helperfun.R")
source("R/preprocessing.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")


# Prediction Part ---------------------------------------------------------

# True FPCs but different cut-off values
m_tru <- list.files(path = paste0(results_wd, "/scen_I_051022/bamlss_tru"))
m_tru_1 <- m_tru[which(as.numeric(substr(m_tru, 2, 4)) < 150)]
m_tru_975 <- m_tru[which(as.numeric(substr(m_tru, 2, 4)) > 149)]
m_tru_95 <- list.files(path = paste0(results_wd, "/scen_I_130922/bamlss_tru"))
m_est_95 <- list.files(path = paste0(results_wd, "/scen_I_130922/bamlss_est"))
m_jmb <- list.files(path = paste0(results_wd, "/scen_I_130922/jmb"))

load(paste0(results_wd, "/scen_I_130922/bamlss_tru/b101.Rdata"))
mfpca_tru <- attr(b_est, "FPCs")

pred_tru_1 <- lapply(m_tru_1, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(results_wd, "/scen_I_051022/bamlss_tru/", x))
  load(paste0(results_wd, "/scen_I_051022/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  vals <- which(mfpca_tru$values > 0)
  newdat <- attach_wfpc(mfpca_tru, d_rirs$data_full, n = length(vals))
  
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
save(pred_tru_1, file = paste0(results_wd, "/pred_tru_1.Rdata"))


pred_tru_975 <- lapply(m_tru_975, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(results_wd, "/scen_I_051022/bamlss_tru/", x))
  load(paste0(results_wd, "/scen_I_051022/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  vals <- which(mfpca_tru$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_tru$values[vals])/sum(mfpca_tru$values[vals]) > 0.975))
  newdat <- attach_wfpc(mfpca_tru, d_rirs$data_full, n = nfpc)
  
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
save(pred_tru_975, file = paste0(results_wd, "/pred_tru_975.Rdata"))


pred_tru_95 <- lapply(m_tru_95, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(results_wd, "/scen_I_130922/bamlss_tru/", x))
  load(paste0(results_wd, "/scen_I_130922/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  vals <- which(mfpca_tru$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_tru$values[vals])/sum(mfpca_tru$values[vals]) > 0.95))
  newdat <- attach_wfpc(mfpca_tru, d_rirs$data_full, n = nfpc)
  
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
save(pred_tru_95, file = paste0(results_wd, "/pred_tru_95.Rdata"))


pred_est_95 <- lapply(m_est_95, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(results_wd, "/scen_I_130922/bamlss_est/", x))
  load(paste0(results_wd, "/scen_I_130922/data/d", substr(x, 2, 4), ".Rdata"))
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Prepare data.frame for the timepoint predictions
  mfpca_est <- attr(b_est, "FPCs")
  vals <- which(mfpca_est$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
  newdat <- na.omit(attach_wfpc(mfpca_est, d_rirs$data_full, n = nfpc))
  
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
save(pred_est_95, file = paste0(results_wd, "/pred_est_95.Rdata"))


pred_jmb <- lapply(m_jmb, function (x) {
  
  # Load the fitted model and the original data
  load(paste0(results_wd, "/scen_I_221019/jmb/", x))
  load(paste0(results_wd, "/scen_I_130922/data/d", substr(x, 5, 7), ".Rdata"))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  
  # Longitudinal fits
  X <- jmb$model_data$X
  Z <- jmb$model_data$Z
  B <- jmb$mcmc$b[[1]]
  mcmc_mu <- do.call(rbind, lapply(1:6, function (dim) {
    tcrossprod(X[[dim]], jmb$mcmc[[paste0("betas", dim)]][[1]]) +
      t(sapply(seq_len(nrow(Z[[dim]])), function (i) {
        Z[[dim]][i, ] %*% B[jmb$model_data$idL[[dim]][i], 
                            (dim - 1)*2 + c(1, 2), ]
      }))
  }))
  
  X_long <- split.data.frame(stats::model.matrix(~obstime*x3,
                                                 data = d_rirs$data_full),
                             d_rirs$data_full$marker)
  Z_long <- split.data.frame(stats::model.matrix(~obstime,
                                                 data = d_rirs$data_full),
                             d_rirs$data_full$marker)
  id_long <- split(d_rirs$data_full$id, d_rirs$data_full$marker)
  mcmc_mu_long <- do.call(rbind, lapply(1:6, function (dim) {
    tcrossprod(X_long[[dim]], jmb$mcmc[[paste0("betas", dim)]][[1]]) +
      t(sapply(seq_len(nrow(Z_long[[dim]])), function (i) {
        Z_long[[dim]][i, ] %*% B[id_long[[dim]][i], 
                                 (dim - 1)*2 + c(1, 2), ]
      }))
  }))
  
  # Survival fits
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                        data = d_rirs$data)[[1]]$knots
  Z <- splineDesign(knots = kn, x = d_rirs$data$obstime, ord = 4, 
                    outer.ok = TRUE)
  X <- jmb$model_data$W_h[unlist(jmb$model_data$idL), , drop = FALSE]
  B <- jmb$mcmc$bs_gammas[[1]]
  Beta <- jmb$mcmc$gammas[[1]]
  mcmc_lambga <- (tcrossprod(Z, B) + tcrossprod(X, Beta))[nodupl_ids, ]
  
  list("lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.025),
                             "Mean" = rowMeans(mcmc_lambga),
                             "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                             probs = 0.975)),
       "alpha" = data.frame("2.5%" = jmb$statistics$CI_low$alphas,
                            "Mean" = jmb$statistics$Mean$alphas,
                            "97.5%" = jmb$statistics$CI_upp$alphas,
                            "marker" = factor(paste0("m", seq_len(6)))),
       "mu" = data.frame("2.5%" = apply(mcmc_mu, 1, quantile, 
                                        probs = 0.025),
                         "Mean" = rowMeans(mcmc_mu),
                         "97.5%" = apply(mcmc_mu, 1, quantile, 
                                         probs = 0.975),
                         "marker" = d_rirs$data$marker),
       "sigma" = data.frame("2.5%" = jmb$statistics$CI_low$sigmas,
                            "Mean" = jmb$statistics$Mean$sigmas,
                            "97.5%" = jmb$statistics$CI_upp$sigmas,
                            "marker" = factor(paste0("m", seq_len(6)))),
       "mu_long" = data.frame("2.5%" = apply(mcmc_mu_long, 1, quantile, 
                                             probs = 0.025),
                              "Mean" = rowMeans(mcmc_mu_long),
                              "97.5%" = apply(mcmc_mu_long, 1, quantile, 
                                              probs = 0.975),
                              "marker" = d_rirs$data_full$marker,
                              "obstime" = d_rirs$data_full$obstime))
  
})
save(pred_jmb, file = paste0(results_wd, "/pred_jmb.Rdata"))



# Simulated Data Part -----------------------------------------------------

d_tru_1 <- paste0("d", substr(m_tru_1, 2, 10))
d_tru_975 <- paste0("d", substr(m_tru_975, 2, 10))
d_tru_95 <- paste0("d", substr(m_tru_95, 2, 10))
d_est_95 <- paste0("d", substr(m_est_95, 2, 10))
d_jmb <- paste0("d", substr(m_jmb, 5, 13))

sdat_tru_1 <- lapply(d_tru_1, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_tru_1, file = paste0(results_wd, "/sdat_tru_1.Rdata"))

sdat_tru_975 <- lapply(d_tru_975, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_tru_975, file = paste0(results_wd, "/sdat_tru_975.Rdata"))

sdat_tru_95 <- lapply(d_tru_95, function (x) {
  
  load(paste0(results_wd, "/scen_I_130922/data/", x))
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  list("lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
       "alpha" = d_rirs$data$alpha[marks],
       "mu" = d_rirs$data[, c("mu", "marker")],
       "sigma" = d_rirs$data$sigma[marks],
       "mu_long" = d_rirs$data_full$mu)
})
save(sdat_tru_95, file = paste0(results_wd, "/sdat_tru_95.Rdata"))

sdat_est_95 <- lapply(d_est_95, function (x) {
  
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
save(sdat_est_95, file = paste0(results_wd, "/sdat_est_95.Rdata"))

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

load(paste0(results_wd, "pred_tru_1.Rdata"))
load(paste0(results_wd, "pred_tru_975.Rdata"))
load(paste0(results_wd, "pred_tru_95.Rdata"))
load(paste0(results_wd, "pred_est_95.Rdata"))
load(paste0(results_wd, "pred_jmb.Rdata"))
load(paste0(results_wd, "sdat_tru_1.Rdata"))
load(paste0(results_wd, "sdat_tru_975.Rdata"))
load(paste0(results_wd, "sdat_tru_95.Rdata"))
load(paste0(results_wd, "sdat_est_95.Rdata"))
load(paste0(results_wd, "sdat_jmb.Rdata"))

r_tru_1 <- sim_results(pred_tru_1, sdat_tru_1, name = "tru_1")
r_tru_975 <- sim_results(pred_tru_975, sdat_tru_975, name = "tru_975")
r_tru_95 <- sim_results(pred_tru_95, sdat_tru_95, name = "tru_95")
r_est_95 <- sim_results(pred_est_95, sdat_est_95, name = "est_95")
r_jmb <- sim_results(pred_jmb, sdat_jmb, name = "jmb")
save(r_tru_1, r_tru_975, r_tru_95, r_est_95, r_jmb, 
     file = paste0(results_wd, "/results.Rdata"))

