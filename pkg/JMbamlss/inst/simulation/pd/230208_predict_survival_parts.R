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
plot_wd <- switch(location,
                  "workstation" = paste0("/run/user/1000/gvfs/smb-share:",
                                         "server=clapton.wiwi.hu-berlin.de,",
                                         "share=volkmana.hub/JMbamlss/",
                                         "simulation/result_plots/scen_II/",
                                         "230208/"),
                  "server_linux" = paste0("~/H:/volkmana.hub/JMbamlss/simul",
                                          "ation/result_plots/scen_II/",
                                          "230208/"))


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


# Define Functions --------------------------------------------------------



# Bamlss survival prediction
sim_bamlss_predict_surv_i <- function(m, wd, model_wd, data_wd, rds = TRUE,
                                      gamma = function(x) 0.3*x, 
                                      lambda = function(t) {
                                        1.65 * t^(0.65) -3
                                      }) {
  
  
  # Load the data set and extract information about it
  if (rds) {
    b_est <- readRDS(paste0(wd, model_wd, m))
    d_rirs <- readRDS(paste0(wd, data_wd, "d", substr(m, 2, 4), ".rds"))
  } else {
    load(paste0(wd, model_wd, m))
    load(paste0(wd, data_wd, "d", substr(m, 2, 4), ".Rdata"))
  }
  nodupl_ids <- which(!duplicated(b_est$model.frame[, c("id", "obstime")]))
  
  
  # Extract MCMC samples to calculate the survival part
  mcmc_lambda <- as.matrix(predict(b_est, model = "lambda", 
                                   FUN = function(x) {x})[nodupl_ids, ])
  mcmc_gamma <- as.matrix(predict(b_est, model="gamma",
                                  FUN = function(x) {x})[nodupl_ids, ])
  mcmc_lambga <- mcmc_gamma + mcmc_lambda
  
  # Hypothetical data to eval gamma coef estimate and lambda on grid
  dat_hypo <- d_rirs$data_hypo[1:101, ]
  dat_hypo$survtime <- seq(0, 1, length.out = 101)
  dat_hypo$x3 <- 1
  
  # Evaluate the Baselinehazard separately
  mcmc_lambda_t <- as.matrix(predict(b_est, dat_hypo, model = "lambda", 
                                     FUN = function(x) {x}))
  mcmc_gamma_int <- as.matrix(predict(b_est, dat_hypo, model = "gamma", 
                                      term = "intercept",
                                      FUN = function(x) {x}))
  mcmc_lambda_tadj <- mcmc_lambda_t + mcmc_gamma_int
  
  
  # Output list
  list("predictions" =  list(
    "lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                         probs = 0.025),
                          "Mean" = rowMeans(mcmc_lambga),
                          "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                          probs = 0.975)),
    "lambda_t" = data.frame("2.5%" = apply(mcmc_lambda_tadj, 1, quantile, 
                                           probs = 0.025),
                            "Mean" = rowMeans(mcmc_lambda_tadj),
                            "97.5%" = apply(mcmc_lambda_tadj, 1, quantile, 
                                            probs = 0.975)),
    "beta_3" = predict(b_est, dat_hypo[1, ], model = "gamma", FUN = c95, 
                       term = 1, intercept = FALSE)
  ),
  "simulations" = list(
    "lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
    "lambda_t" = lambda(dat_hypo$survtime), 
    "beta_3" = gamma(dat_hypo[1, "x3"])
  ))
  
}

#' Simulation Helper Function - Predict the Results for bamlss-Models
#' 
#' This function takes all the models listed in a folder and predicts the fit.
#' 
#' @param m Vector containing all the file names of the models. 
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_bamlss_predict_surv <- Vectorize(sim_bamlss_predict_surv_i, 
                                     vectorize.args = "m", SIMPLIFY = FALSE)


#' Simulation Helper Function - Evaluate the Simulation for JMbamlss Setting
#' 
#' This function evaluates the results for a given folder of JMbamlss model 
#' fits.
#' 
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param name Name for description of the simulation setting.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_jmbamlss_surv <- function(wd, model_wd, data_wd, name, rds = TRUE) {
  
  models <- list.files(path = paste0(wd, model_wd))
  if (rds) {
    models <- models[grep("\\.rds", models)]
  } else {
    rmv <- grep("\\.rds", models)
    if(length(rmv) > 0) {
      models <- models[-grep("\\.rds", models)]
    }
  }
  list_to_compare <- sim_bamlss_predict_surv(models, wd, model_wd, data_wd, rds)
  
}


sim_surv_results_list <- function(result_list, name) {
  
  lapply(result_list, function(res) {
      
    est <- res$predictions
    sim <- res$simulations
    
    # Bias, MSE, Coverage
    eval_lambga <- data.frame(
      type = c("Bias", "MSE", "Coverage"),
      model = name,
      predictor = "lambga",
      t = "all",
      value = c(mean(-(est$lambga$Mean - sim$lambga)),
                mean((est$lambga$Mean - sim$lambga)^2),
                mean(est$lambga[, 1] < sim$lambga & 
                       est$lambga[, 3] > sim$lambga)))
    
    eval_beta_3 = data.frame(
      type = rep(c("Bias", "MSE", "Coverage")),
      model = name,
      predictor = "beta_3",
      t = "all",
      value = c(-(est$beta_3$Mean - sim$beta_3),
                (est$beta_3$Mean - sim$beta_3)^2,
                as.numeric(est$beta_3[, 1] < sim$beta_3 & 
                             est$beta_3[, 3] > sim$beta_3)))
    
    eval_lambda_t <- data.frame(
      type = rep(c("Bias", "MSE", "Coverage"), each = 101),
      model =  name,
      predictor = "lambda_t",
      t = rep(seq(0, 1, by = 0.01), times = 3),
      value = c(sim$lambda_t - est$lambda_t$Mean,
                (sim$lambda_t - est$lambda_t$Mean)^2,
                as.numeric(est$lambda_t[, 1] < sim$lambda_t &
                             est$lambda_t[, 3] > sim$lambda_t)))
    
    
    
    
    rbind(eval_lambga, eval_beta_3, eval_lambda_t)
    
  })
}
sim_surv_results <- function(result_list, name) {
  it_list <- sim_surv_results_list(result_list, name)
  do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)), it_list))
}


sim_jmb_predict_surv_i <- function(m, wd, model_wd, data_wd, 
                                   gamma = function(x) 0.3*x, 
                                   lambda = function(t) {
                                     1.65 * t^(0.65) -3
                                   }, rds = TRUE) {
  
  # Load the fitted model and the original data
  if (rds) {
    jmb <- readRDS(paste0(wd, model_wd, m))
    d_rirs <- readRDS(paste0(wd, data_wd, "d", substr(m, 4, 6), ".rds"))
  } else {
    load(paste0(wd, model_wd, m))
    load(paste0(wd, data_wd, "d", substr(m, 4, 6), ".Rdata"))
  }
  
  nodupl_ids <- which(!duplicated(d_rirs$data[, c("id", "obstime")]))
  n_dim <- length(levels(d_rirs$data$marker))
  marks <- which(!duplicated(d_rirs$data$marker))
  
  
  # Survival fits
  kn <- mgcv::smoothCon(mgcv::s(survtime, k = 20, bs = "ps"), 
                        data = d_rirs$data)[[1]]$knots
  Z <- splineDesign(knots = kn, x = d_rirs$data$obstime, ord = 4, 
                    outer.ok = TRUE)
  X <- jmb$model_data$W_h[unlist(jmb$model_data$idL), , drop = FALSE]
  B <- jmb$mcmc$bs_gammas[[1]]
  Beta <- jmb$mcmc$gammas[[1]]
  mcmc_lambga <- (tcrossprod(Z, B) + tcrossprod(X, Beta))[nodupl_ids, ]
  
  # Hypothetical data to eval gamma coef estimate and lambda on grid
  dat_hypo <- d_rirs$data_hypo[1:101, ]
  dat_hypo$survtime <- seq(0, 1, length.out = 101)
  
  # Lambda fit
  Z_t <- splineDesign(knots = kn, x = dat_hypo$survtime, ord = 4, 
                      outer.ok = TRUE)
  mcmc_lambda_t <- tcrossprod(Z_t, B)
  
  list("predictions" = list(
    "lambga" = data.frame("2.5%" = apply(mcmc_lambga, 1, quantile, 
                                         probs = 0.025),
                          "Mean" = rowMeans(mcmc_lambga),
                          "97.5%" = apply(mcmc_lambga, 1, quantile, 
                                          probs = 0.975)),
    "lambda_t" = data.frame("2.5%" = apply(mcmc_lambda_t, 1, quantile, 
                                           probs = 0.025),
                            "Mean" = rowMeans(mcmc_lambda_t),
                            "97.5%" = apply(mcmc_lambda_t, 1, quantile, 
                                            probs = 0.975)),
    "beta_3" = data.frame("2.5%" = jmb$statistics$CI_low$gammas,
                          "Mean" = jmb$statistics$Mean$gammas,
                          "97.5%" = jmb$statistics$CI_upp$gammas)),
    "simulations" = list(
      "lambga" = rowSums(d_rirs$data[nodupl_ids, c("lambda", "gamma")]),
      "lambda_t" = lambda(dat_hypo$survtime), 
      "beta_3" = gamma(1)
    ))
  
  
}

#' Simulation Helper Function - Predict the Results for bamlss-Models
#' 
#' This function takes all the models listed in a folder and predicts the fit.
#' 
#' @param m Vector containing all the file names of the models. 
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_jmb_predict_surv <- Vectorize(sim_jmb_predict_surv_i, 
                                  vectorize.args = "m", SIMPLIFY = FALSE)

#' Simulation Helper Function - Evaluate the Simulation for JMbamlss Setting
#' 
#' This function evaluates the results for a given folder of JMbamlss model 
#' fits.
#' 
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param name Name for description of the simulation setting.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_jmb_surv <- function(wd, model_wd, data_wd, name, rds = TRUE) {
  
  models <- list.files(path = paste0(wd, model_wd))
  if (rds) {
    models <- models[grep("\\.rds", models)]
  } else {
    rmv <- grep("\\.rds", models)
    if(length(rmv) > 0) {
      models <- models[-grep("\\.rds", models)]
    }
  }
  list_to_compare <- sim_jmb_predict_surv(models, wd, model_wd, data_wd, rds)
  
}

# Apply Functions ---------------------------------------------------------


scenII_A2 <- sim_jmbamlss_surv(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "A2/", data_wd = "data/", name = "A", rds = TRUE)
saveRDS(scenII_A2, file = paste0(server_wd, "scen_II_230117/surv_pred_A2.rds"))

scenII_E <- sim_jmbamlss_surv(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "E/", data_wd = "data/", name = "E", rds = TRUE)
saveRDS(scenII_E, file = paste0(server_wd, "scen_II_230117/surv_pred_E.rds"))

scenII_F2 <- sim_jmbamlss_surv(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "F2/", data_wd = "data/", name = "F", rds = TRUE)
saveRDS(scenII_F2, file = paste0(server_wd, "scen_II_230117/surv_pred_F2.rds"))

scenII_L2 <- sim_jmbamlss_surv(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "L2/", data_wd = "data/", name = "L", rds = TRUE)
saveRDS(scenII_L2, file = paste0(server_wd, "scen_II_230117/surv_pred_L2.rds"))

scenII_JMB <- sim_jmb_surv(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"), 
  model_wd = "JMB/", data_wd = "data/", name = "JMB", rds = TRUE)
saveRDS(scenII_JMB, file = paste0(server_wd,
                                  "scen_II_230117/surv_pred_JMB.rds"))




# Look at Evaluations -----------------------------------------------------

scenII_A2 <- readRDS(paste0(server_wd, "scen_II_230117/surv_pred_A2.rds"))
res_A <- sim_surv_results(scenII_A2, "A")

scenII_F2 <- readRDS(paste0(server_wd, "scen_II_230117/surv_pred_F2.rds"))
res_F <- sim_surv_results(scenII_F2, "F")

scenII_L2 <- readRDS(paste0(server_wd, "scen_II_230117/surv_pred_L2.rds"))
res_L <- sim_surv_results(scenII_L2, "L")

scenII_E <- readRDS(paste0(server_wd, "scen_II_230117/surv_pred_E.rds"))
res_E <- sim_surv_results(scenII_E, "E")

scenII_JMB <- readRDS(paste0(server_wd, "scen_II_230117/surv_pred_JMB.rds"))
res_JMB <- sim_surv_results(scenII_JMB, "JMB")


eval_dat <- rbind(res_A, res_E, res_F, res_L, res_JMB) %>%
  mutate(model = factor(model, levels = c("E", "L", "F", "A", "JMB"),
                        labels = c("True", "Uncens", "Cens", "Cens95", "JMB")))

# Same Graphs can be recreated
ggplot(eval_dat %>% filter(predictor == "lambga", type == "MSE"),
       aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("MSE: Lambda + Gamma")

# Gamma coefficient estimate
ggplot(eval_dat %>% filter(predictor == "beta_3", type == "MSE"),
       aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("MSE: Gamma (Coef)")
ggsave(paste0(plot_wd, "gamma_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(eval_dat %>% filter(predictor == "beta_3", type == "Bias"),
       aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. Bias") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("Bias: Gamma (Coef)")
ggsave(paste0(plot_wd, "gamma_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

eval_dat %>% 
  filter(predictor == "beta_3", type == "Coverage") %>%
  group_by(model) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)



# Lambda evaluated over t
ggplot(eval_dat %>% filter(predictor == "lambda_t", type == "MSE", 
                      t %in% seq(0, 1, by = 0.1)),
       aes(x = t, y = value, color = model)) +
  geom_boxplot()

ggplot(data = eval_dat %>% 
         filter(predictor == "lambda_t",
                type == "MSE",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Lambda (by t)")+
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "lambdat_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")
# Cutoff
ggplot(data = eval_dat %>% 
         filter(predictor == "lambda_t",
                type == "MSE",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Lambda (by t)")+
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw()
ggsave(paste0(plot_wd, "lambdat_mse_cutoff.pdf"), device = "pdf", width = 8,
       height = 6, units = "in")

ggplot(data = eval_dat %>% 
         filter(predictor == "lambda_t",
                type == "Bias",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  ggtitle("Bias: Lambda (by t)")+
  ylab("Emp. Bias") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "lambdat_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% 
         filter(predictor == "lambda_t",
                type == "Coverage") %>%
         group_by(model, t) %>%
         summarise(Cov = mean(value), .groups = "keep") %>%
         ungroup(),
       aes(y = Cov, x = as.numeric(t), col = model)) +
  geom_line() +
  ggtitle("Coverage: Lambda (by t)")+
  ylab("Emp. Coverage") +
  xlab("t") +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "lambdat_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")



# Plot the estimated Baseline Hazard --------------------------------------

dat_plot <- rbind(
  do.call(rbind, 
          Map(cbind, it = seq_along(scenII_A2),
              value = lapply(scenII_A2, 
                             function(x) x$predictions$lambda_t$Mean))),
  do.call(rbind, 
          Map(cbind, it = seq_along(scenII_E),  
              value = lapply(scenII_E, 
                             function(x) x$predictions$lambda_t$Mean))),
  do.call(rbind, 
          Map(cbind, it = seq_along(scenII_F2), 
              value = lapply(scenII_F2,
                             function(x) x$predictions$lambda_t$Mean))),
  do.call(rbind, 
          Map(cbind, it = seq_along(scenII_L2), 
              value = lapply(scenII_L2, 
                             function(x) x$predictions$lambda_t$Mean))),
  do.call(rbind, 
          Map(cbind, it = seq_along(scenII_JMB), 
              value = lapply(scenII_JMB, 
                             function(x) x$predictions$lambda_t$Mean)))
) %>%
  as.data.frame() %>%
  mutate(model = factor(rep(c("A", "E", "F", "L", "JMB"), each = 101*100)),
         obstime = rep(seq(0, 1, by = 0.01), 5*100),
         type = factor("Est")) %>%
  rbind(data.frame(it = 0, 
                   value = rep(scenII_E$b100.rds$simulations$lambda_t, 5),
                   model = rep(c("A", "E", "F", "L", "JMB"), each = 101),
                   obstime = seq(0, 1, by = 0.01), type = factor("Tru"))) %>%
  mutate(model = factor(model, levels = c("E", "L", "F", "A", "JMB"),
                        labels = c("True", "Uncens", "Cens", "Cens95", "JMB")))

ggplot(dat_plot,
       aes(x = obstime, y = value, group = it, color = type, alpha = type)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~model) +
  ggtitle("Estimated Baseline Hazards") +
  labs(color = NULL, alpha = NULL, x = "Time", y = "Baseline Hazards")
ggsave(paste0(plot_wd, "baselineH.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")
