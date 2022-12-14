library(bamlss)
library(JMbayes2)
library(funData)
library(tidyverse)
source("R/preprocessing.R")
source("R/eval_mfun.R")

# Extract Simulation Results ----------------------------------------------

wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation/scen_I_130922")
m_est <- list.files(path = paste0(wd, "/bamlss_est"))
m_tru <- list.files(path = paste0(wd, "/bamlss_tru"))
m_jmb <- list.files(path = paste0(wd, "/jmb"))
l_dat <- list.files(path = paste0(wd, "/data"))


full <- substr(m_est, 2, 4)
m_tru <- m_tru[100:199 %in% full]
m_jmb <- m_jmb[100:199 %in% full]
l_dat <- l_dat[100:199 %in% full]

res_est <- lapply(m_est, function (x) {
  load(paste0(wd, "/bamlss_est/", x))
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

res_tru <- lapply(m_tru, function (x) {
  
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
  mfpca_est <- attr(b_est, "FPCs")
  vals <- which(mfpca_est$values > 0)
  nfpc <- min(which(
    cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
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

res_jmb <- lapply(m_jmb, function (x) {
  
  # Load the fitted model and the original data
  load(paste0(wd, "/jmb/", x))
  load(paste0(wd, "/data/d", substr(x, 5, 7), ".Rdata"))
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
   
  # Does not work
  # pred <- predict(jmb, newdata = list(newdataL = d_rirs$data,
  #                             newdataE = d_rirs$data_short[1:150, ]))
  
})

sim_dat <- lapply(l_dat, function (x) {
  
  load(paste0(wd, "/data/", x))
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

# Alpha, sigma, mu sollte man getrennt nach Markern betrachten
results <- mapply(function (est, tru, jmb, sim) {

  # Bias, MSE, Coverage
  eval_lambga <- data.frame(
    type = rep(c("Bias", "MSE", "Coverage"), each = 3),
    model = c("est","tru", "jmb") ,
    predictor = "lambga",
    marker = "all",
    t = "all",
    value = c(mean(est$lambga$Mean - sim$lambga),
              mean(tru$lambga$Mean - sim$lambga),
              mean(jmb$lambga$Mean - sim$lambga),
              mean((est$lambga$Mean - sim$lambga)^2),
              mean((tru$lambga$Mean - sim$lambga)^2),
              mean((jmb$lambga$Mean - sim$lambga)^2),
              mean(est$lambga[, 1] < sim$lambga & est$lambga[, 3] > sim$lambga),
              mean(tru$lambga[, 1] < sim$lambga & tru$lambga[, 3] > sim$lambga),
              mean(jmb$lambga[, 1] < sim$lambga & jmb$lambga[, 3] > sim$lambga)))
  
  eval_alpha = data.frame(
    type = rep(c("Bias", "MSE", "Coverage"), each = 3*6),
    model = rep(c("est","tru", "jmb"), each = 6),
    predictor = "alpha",
    marker = paste0("m", 1:6),
    t = "all",
    value = c(est$alpha$Mean - sim$alpha,
              tru$alpha$Mean - sim$alpha,
              jmb$alpha$Mean - sim$alpha,
              (est$alpha$Mean - sim$alpha)^2,
              (tru$alpha$Mean - sim$alpha)^2,
              (jmb$alpha$Mean - sim$alpha)^2,
              as.numeric(est$alpha[, 1] < sim$alpha & 
                           est$alpha[, 3] > sim$alpha),
              as.numeric(tru$alpha[, 1] < sim$alpha &
                           tru$alpha[, 3] > sim$alpha),
              as.numeric(jmb$alpha[, 1] < sim$alpha & 
                           jmb$alpha[, 3] > sim$alpha)))
  
  eval_sigma = data.frame(
    type = rep(c("Bias", "MSE", "Coverage"), each = 3*6),
    model = rep(c("est","tru", "jmb"), each = 6),
    predictor = "sigma",
    marker = paste0("m", 1:6),
    t = "all",
    value = c(est$sigma$Mean - sim$sigma,
              tru$sigma$Mean - sim$sigma,
              jmb$sigma$Mean - sim$sigma,
              (est$sigma$Mean - sim$sigma)^2,
              (tru$sigma$Mean - sim$sigma)^2,
              (jmb$sigma$Mean - sim$sigma)^2,
              as.numeric(est$sigma[, 1] < sim$sigma & 
                           est$sigma[, 3] > sim$sigma),
              as.numeric(tru$sigma[, 1] < sim$sigma &
                           tru$sigma[, 3] > sim$sigma),
              as.numeric(jmb$sigma[, 1] < sim$sigma & 
                           jmb$sigma[, 3] > sim$sigma)))
  #browser()
  eval_mu <- data.frame(
    type = rep(c("Bias", "MSE", "Coverage"), each = 3*6),
    model =  rep(c("est","tru", "jmb"), each = 6),
    predictor = "mu",
    marker = paste0("m", 1:6),
    t = "all",
    value = c(mapply(function (e, s) {
                  mean(e - s)
                }, e = split(est$mu$Mean, est$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (t, s) {
                  mean(t - s)
                }, t = split(tru$mu$Mean, tru$mu$marker), 
                s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (j, s) {
                mean(j - s)
              }, j = split(jmb$mu$Mean, jmb$mu$marker), 
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (e, s) {
                mean((e - s)^2)
              }, e = split(est$mu$Mean, est$mu$marker), 
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (t, s) {
                mean((t - s)^2)
              }, t = split(tru$mu$Mean, tru$mu$marker), 
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (j, s) {
                mean((j - s)^2)
              }, j = split(jmb$mu$Mean, jmb$mu$marker), 
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (l, u, s) {
                mean(as.numeric(l < s & u > s))
              }, l = split(est$mu[, 1], est$mu$marker), 
              u = split(est$mu[, 3], est$mu$marker),
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (l, u, s) {
                mean(as.numeric(l < s & u > s))
              }, l = split(tru$mu[, 1], tru$mu$marker), 
              u = split(tru$mu[, 3], tru$mu$marker),
              s = split(sim$mu$mu, sim$mu$marker)),
              mapply(function (l, u, s) {
                mean(as.numeric(l < s & u > s))
              }, l = split(jmb$mu[, 1], jmb$mu$marker), 
              u = split(jmb$mu[, 3], jmb$mu$marker),
              s = split(sim$mu$mu, sim$mu$marker))))
  
  sim_marker <- split(sim$mu_long, est$mu_long$marker)
  eval_mu_long <- data.frame(
    type = rep(c("Bias", "MSE", "Coverage"), each = 3*6*101),
    model =  rep(c("est", "tru", "jmb"), each = 6*101),
    predictor = "mu_long",
    marker = paste0("m", 1:6),
    t = rep(rep(seq(0, 1, by = 0.01), each = 6), times = 3),
    value = c(c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e$Mean[same_t] - s[same_t])
      }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e$Mean[same_t] - s[same_t])
      }, e = split(tru$mu_long, tru$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e$Mean[same_t] - s[same_t])
      }, e = split(jmb$mu_long, jmb$mu_long$marker), s = sim_marker)
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
        mean((e$Mean[same_t] - s[same_t])^2)
      }, e = split(tru$mu_long, tru$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean((e$Mean[same_t] - s[same_t])^2)
      }, e = split(jmb$mu_long, jmb$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
      }, e = split(est$mu_long, est$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
      }, e = split(tru$mu_long, tru$mu_long$marker), s = sim_marker)
    })),
    c(sapply(seq(0, 1, by = 0.01), function (t) {
      mapply(function(e, s) {
        same_t <- e$obstime == t
        mean(e[same_t, 1] < s[same_t] & e[same_t, 3] > s[same_t])
      }, e = split(jmb$mu_long, jmb$mu_long$marker), s = sim_marker)
    }))))
  
  rbind(eval_lambga, eval_alpha, eval_sigma, eval_mu, eval_mu_long)

}, est = res_est, tru = res_tru, jmb = res_jmb, sim = sim_dat, SIMPLIFY = FALSE)

ggplot(data = do.call(rbind, results) %>% filter(predictor == "lambga",
                                                 type == "Bias"),
       aes(y = value, x = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5))
ggplot(data = do.call(rbind, results) %>% filter(predictor == "lambga", 
                                                 type == "MSE"),
       aes(y = value, x = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5))
ggplot(data = do.call(rbind, results) %>% filter(predictor == "lambga", 
                                                 type == "Coverage" ),
       aes(y = value, x = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "alpha",
                                                 type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5))
ggplot(data = do.call(rbind, results) %>% filter(predictor == "alpha",
                                                 type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5))
ggplot(data = do.call(rbind, results) %>% filter(predictor == "alpha",
                                                 type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu",
                                                 type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu",
                                                 type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu",
                                                 type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "sigma",
                                                 type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "sigma",
                                                 type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "sigma",
                                                 type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot()
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu_long",
                                                 type == "Bias",
                                                 marker %in% c("m1", "m2", "m3"),
                                                 t %in% seq(0, 1, by = 0.02)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(-0.25, 0.25))
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu_long",
                                                 type == "Bias",
                                                 marker %in% c("m1", "m2", "m3"),
                                                 t %in% seq(0, 1, by = 0.02),
                                                 model != "est"),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free")
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu_long",
                                                 type == "MSE"),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free")
ggplot(data = do.call(rbind, results) %>% filter(predictor == "mu_long",
                                                 type == "Coverage",
                                                 marker %in% c("m1", "m2", "m3"),
                                                 t %in% seq(0, 1, by = 0.02),
                                                 model != "est"),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free")


