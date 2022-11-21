sim_results <- function(result_list, dat_list, name = "est") {
  
  mapply(function (est, sim) {
    
    # Bias, MSE, Coverage
    eval_lambga <- data.frame(
      type = c("Bias", "MSE", "Coverage"),
      model = name,
      predictor = "lambga",
      marker = "all",
      t = "all",
      value = c(mean(-(est$lambga$Mean - sim$lambga)),
                mean((est$lambga$Mean - sim$lambga)^2),
                mean(est$lambga[, 1] < sim$lambga & 
                       est$lambga[, 3] > sim$lambga)))
    
    eval_alpha = data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model = name,
      predictor = "alpha",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(-(est$alpha$Mean - sim$alpha),
                (est$alpha$Mean - sim$alpha)^2,
                as.numeric(est$alpha[, 1] < sim$alpha & 
                             est$alpha[, 3] > sim$alpha)))
    
    eval_sigma = data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model = name,
      predictor = "sigma",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(-(est$sigma$Mean - sim$sigma),
                (est$sigma$Mean - sim$sigma)^2,
                as.numeric(est$sigma[, 1] < sim$sigma & 
                             est$sigma[, 3] > sim$sigma)))
    
    eval_mu <- data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6),
      model =  name,
      predictor = "mu",
      marker = paste0("m", 1:6),
      t = "all",
      value = c(mapply(function (e, s) {
                  mean(s - e)
                }, e = split(est$mu$Mean, est$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (e, s) {
                  mean((e - s)^2)
                }, e = split(est$mu$Mean, est$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
                mapply(function (l, u, s) {
                  mean(as.numeric(l < s & u > s))
                }, l = split(est$mu[, 1], est$mu$marker), 
                u = split(est$mu[, 3], est$mu$marker),
                s = split(sim$mu$mu, sim$mu$marker))))
    
    sim_marker <- split(sim$mu_long, est$mu_long$marker)
    eval_mu_long <- data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 6*101),
      model =  name,
      predictor = "mu_long",
      marker = paste0("m", 1:6),
      t = rep(rep(seq(0, 1, by = 0.01), each = 6), times = 3),
      value = c(c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(-(e$Mean[same_t] - s[same_t]))
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean((e$Mean[same_t] - s[same_t])^2)
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                })),
                c(sapply(seq(0, 1, by = 0.01), function (t) {
                  mapply(function(e, s) {
                    same_t <- e$obstime == t
                    mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
                  }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
                }))))
    
    rbind(eval_lambga, eval_alpha, eval_sigma, eval_mu, eval_mu_long)
    
  }, est = result_list, sim = dat_list, SIMPLIFY = FALSE)
  
}


sim_bamlss_predict_i <- function(i, results_wd, m_setting = "/scen_I_130922/",
                                 d_setting = "/scen_I_130922/",
                                 folder = "bamlss_est/") {
  
  # Load the data set and extract information about it
  load(paste0(results_wd, m_setting, folder, i))
  load(paste0(results_wd, d_setting, "data/d", substr(i, 2, 4), ".Rdata"))
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
}

#' Simulation Helper Function - Predict the Results for bamlss-Models
#' 
#' This function takes all the models listed in a folder and predicts the fit.
#' 
#' @param m Vector containing all the file names of the models. 
#' @param results_wd Path to simulations folder.
#' @param m_setting Simulation setting folder where the models are saved.
#' @param d_setting Simulation data folder.
#' @param folder Model folder where the models are saved (for different
#'   specifications).
sim_bamlss_predict <- Vectorize(sim_bamlss_predict_i, vectorize.args = "i", 
                                SIMPLIFY = FALSE)


sim_jmb_predict_i <- function(i, results_wd, m_setting = "/scen_I_221019/",
                              d_setting = "/scen_I_130922/",
                              folder = "jmb/") {
  
  # Load the fitted model and the original data
  load(paste0(results_wd, m_setting, folder, i))
  load(paste0(results_wd, d_setting, "data/d", substr(i, 5, 7), ".Rdata"))
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
  
}

#' Simulation Helper Function - Predict the Results for JMbayes-Models
#' 
#' This function takes all the models listed in a folder and predicts the fit.
#' 
#' @param m Vector containing all the file names of the models. 
#' @param results_wd Path to simulations folder.
#' @param m_setting Simulation setting folder where the models are saved.
#' @param d_setting Simulation data folder.
#' @param folder Model folder where the models are saved (for different
#'   specifications).
sim_jmb_predict <- Vectorize(sim_jmb_predict_i, vectorize.args = "i", 
                                SIMPLIFY = FALSE)
