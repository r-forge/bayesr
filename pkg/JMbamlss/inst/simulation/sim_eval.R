
# Extract Simulation Results ----------------------------------------------

wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation/scen_I_130922")
m_est <- list.files(path = paste0(wd, "/bamlss_est"))
m_tru <- list.files(path = paste0(wd, "/bamlss_tru"))
m_jmb <- list.files(path = paste0(wd, "/jmb"))

res_est <- lapply(m_est, function (x) {
  load(paste0(wd, "/bamlss_est/", x))
  ids <- which(!duplicated(b_est$model.frame$id))
  id_marks <- which(!duplicated(b_est$model.frame[, c("id", "marker")]))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  list("lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.025),
                             "Mean" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.5),
                             "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                             probs = 0.975)),
       "alpha" = cbind(predict(b_est, model = "alpha", FUN = c95)[id_marks, ],
                       data.frame("marker" = b_est$model.frame$marker[id_marks])),
       "mu" = cbind(predict(b_est, model = "mu", FUN = c95), 
                    data.frame("marker" = b_est$model.frame$marker)),
       "sigma" = cbind(predict(b_est, model = "sigma", FUN = c95)[id_marks,],
                       data.frame("marker" = b_est$model.frame$marker)))
})

res_tru <- lapply(m_tru, function (x) {
  
  # Load the data set and extract information about it
  load(paste0(wd, "/bamlss_tru/", x))
  ids <- which(!duplicated(b_est$model.frame$id))
  id_marks <- which(!duplicated(b_est$model.frame[, c("id", "marker")]))
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  list("lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.025),
                             "Mean" = apply(mcmc_lambga, 1, quantile, 
                                            probs = 0.5),
                             "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                             probs = 0.975)),
       "alpha" = cbind(predict(b_est, model = "alpha", FUN = c95)[id_marks, ],
                       data.frame("marker" = b_est$model.frame$marker[id_marks])),
       "mu" = cbind(predict(b_est, model = "mu", FUN = c95), 
                    data.frame("marker" = b_est$model.frame$marker)),
       "sigma" = cbind(predict(b_est, model = "sigma", FUN = c95)[id_marks,],
                       data.frame("marker" = b_est$model.frame$marker)))
})

res_jmb <- lapply(m_jmb, function (x) {
  
  # Load the fitted model and the original data
  load(paste0(wd, "/jmb/", x))
  load(paste0(wd, "/data/d", substr(x, 5, 7), ".Rdata"))
  predict(jmb, newdata = )
  
})