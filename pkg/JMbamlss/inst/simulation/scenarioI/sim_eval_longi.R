library(bamlss)
library(tidyverse)
library(manipulate)

predict_sim_longi <- function(scenario, iter) {
  load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
              "share=volkmana.hub/JMbamlss/simulation/", scenario, 
              "/bamlss_tru/b", iter, ".Rdata"))
  b_tru <- b_est
  load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
              "share=volkmana.hub/JMbamlss/simulation/", scenario, 
              "/bamlss_est/b", iter, ".Rdata"))
  load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
              "share=volkmana.hub/JMbamlss/simulation/", scenario,
              "/jmb/jmb_", iter, ".Rdata"))
  load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
              "share=volkmana.hub/JMbamlss/simulation/", scenario,
              "/data/d", iter, ".Rdata"))
  
  newdat <- d_rirs$data
  newdat$yest <- predict(b_est, newdata = newdat, model = "mu")
  newdat$ytru <- predict(b_tru, newdata = newdat, model = "mu")
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
  newdat$yjmb <- rowMeans(mcmc_mu)
  
  newdat
  
}

fit_short <- predict_sim_longi(scenario = "scen_I_130922", iter = 120)
manipulate(
  ggplot(fit_short %>% filter(id == subj.new), aes(x = obstime, y = yest)) +
    geom_line(col = "red") +
    facet_grid(~marker) +
    geom_point(aes(y =y)) +
    geom_line(aes(y = yjmb), col = "blue") +
    geom_line(aes(y = ytru), col = "green"),
  subj.new = slider(1, 150)
)
fit_long <- sim_longi(scenario = "scen_I", iter = 120)
manipulate(
  ggplot(fit_long %>% filter(id == subj.new), aes(x = obstime, y = yest)) +
    geom_line(col = "red") +
    facet_grid(~marker) +
    geom_point(aes(y =y)) +
    geom_line(aes(y = yjmb), col = "blue") +
    geom_line(aes(y = ytru), col = "green"),
  subj.new = slider(1, 150)
)


# Why are the models not flexible enough? ---------------------------------



load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_tru",
            "/b101.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/jmb/jmb_10",
            "1.Rdata"))
newdat <- b_est$model.frame  %>%
  mutate(y = b_est$model.frame[, 1][, 3]) %>%
  filter(id == 19)
newdat$yhat <- predict(b_est, newdata = newdat, model = "mu")
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
newdat$yjmb <- rowMeans(mcmc_mu)[unlist(jmb$model_data$idL) == 19]

ggplot(newdat, aes(x = obstime, y = yhat)) +
  geom_line(col = "red") +
  facet_grid(~marker) +
  geom_point(aes(y =y)) +
  geom_line(aes(y = yjmb), col = "blue")

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/bamlss_tru",
            "/b101.Rdata"))
newdat$yhat2 <- predict(b_est, newdata = newdat, model = "mu")
ggplot(newdat, aes(x = obstime, y = yhat)) +
  geom_line(col = "red") +
  facet_grid(~marker) +
  geom_point(aes(y =y)) +
  geom_line(aes(y = yjmb), col = "blue") +
  geom_line(aes(y = yhat2), col = "green")


# Here is the comparison of all 12 FPCs with JMB
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/bamlss_tru",
            "/b100.Rdata"))
newdat <- b_est$model.frame  %>%
  mutate(y = b_est$model.frame[, 1][, 3])
newdat$yhat <- predict(b_est, newdata = newdat, model = "mu")
newdat$yhat2 <- predict(b_est, newdata = newdat, model = "mu", 
                        what = "parameters")
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/jmb/",
            "jmb_100.Rdata"))
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
newdat$yjmb <- rowMeans(mcmc_mu)


manipulate(
  ggplot(newdat %>% filter(id == subj.new), aes(x = obstime, y = yhat)) +
    geom_line(col = "red") +
    facet_grid(~marker) +
    geom_point(aes(y =y)) +
    geom_line(aes(y = yhat2), col = "green") +
    geom_line(aes(y = yjmb), col = "blue"),
  subj.new = slider(1, 150)
)
