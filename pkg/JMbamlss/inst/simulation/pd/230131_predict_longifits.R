
# Compare different fits --------------------------------------------------


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
                                         "simulation/result_plots/scen_II/"),
                  "server_linux" = paste0("~/H:/volkmana.hub/JMbamlss/simul",
                                          "ation/result_plots/scen_II/"))


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

set.seed(1520)
ids <- sample(1:300, size = 5)

# New functions for predicting specific observations ----------------------

sim_predict_longobs_i <- function(obs, m, wd, model_wd, data_wd, rds = TRUE) {
  
  
  # Load the data set and extract information about it
  if (rds) {
    b_est <- readRDS(paste0(wd, model_wd, m))
    d_rirs <- readRDS(paste0(wd, data_wd, "d", substr(m, 2, 4), ".rds"))
  } else {
    load(paste0(wd, model_wd, m))
    load(paste0(wd, data_wd, "d", substr(m, 2, 4), ".Rdata"))
  }
  marks <- which(!duplicated(b_est$model.frame$marker))
  
  # Predict only specific observations
  keep <- b_est$model.frame$id %in% obs
  d_rirs$data_full <- d_rirs$data_full %>% filter(id %in% obs)
  
  # Output list
  list("predictions" =  list(
    "mu" = cbind(predict(b_est, model = "mu", FUN = c95)[keep, ], 
                 data.frame("id" = b_est$model.frame$id[keep],
                            "marker" = b_est$model.frame$marker[keep],
                            "obstime" = b_est$model.frame$obstime[keep])),
    "mu_long" = cbind(predict(b_est, model = "mu", FUN = c95, 
                              newdata = d_rirs$data_full),
                      data.frame("id" = d_rirs$data_full$id,
                                 "marker" = d_rirs$data_full$marker,
                                 "obstime" = d_rirs$data_full$obstime))
  ),
  "simulations" = list(
    "mu" = d_rirs$data[keep, c("id", "marker", "obstime", "mu", "y")],
    "mu_long" = d_rirs$data_full[, c("id", "marker", "obstime", "mu")]
  ))
  
}

#' Simulation Helper Function - Predict the Longitudinal Fits
#' 
#' This function takes all the models listed in a folder and predicts the 
#' longitudinal fit for specific observations.
#' 
#' @param obs Vector containing the observations for which longitudinal 
#'   predictions are to be made.
#' @param m Vector containing all the file names of the models. 
#' @param wd Path to simulations folder.
#' @param model_wd Simulation setting folder where the models are saved.
#' @param data_wd Simulation data folder.
#' @param rds Objects are saved as .rds files (for backwards compatibility when
#'   .Rdata files were used). Defaults to TRUE.
sim_predict_longobs <- Vectorize(sim_predict_longobs_i, vectorize.args = "m", 
                                 SIMPLIFY = FALSE)




# Different Fits with Different Seeds -------------------------------------


fit_F <- sim_predict_longobs(obs = ids, paste0("b", 100:131, ".Rdata"), 
                             wd = paste0(server_wd, "scen_II_221209/"),
                             model_wd = "F/", data_wd = "data/", rds = FALSE)

dat_F <- do.call(rbind, Map(cbind, it = seq_along(fit_F), 
                            lapply(fit_F, function(x) {
  x$predictions$mu_long[, c("id", "marker", "obstime", "Mean")]
})))

dat_F <- rbind(dat_F, fit_F[[1]]$simulations$mu_long %>%
                 mutate(Mean = mu, mu = NULL, it = 0)) %>%
  mutate(src = factor(ifelse(it == 0, "TRUE", "EST")))


ggplot(dat_F, aes(x = obstime, y = Mean, group = it, color = src, 
                  alpha = src)) + 
  geom_line() +
  facet_grid(id ~ marker) +
  scale_alpha_manual(values = c(0.6, 1)) +
  theme_bw() +
  labs(y = expression(mu(t)~","~hat(mu)(t)), alpha = NULL, color = NULL) +
  ggtitle("Different Fits based on Different Seeds")



# Different Fits with Different Seeds and Different FPCs ------------------

fit_L <- sim_predict_longobs(obs = ids, paste0("b", 100:199, ".Rdata"), 
                             wd = paste0(server_wd, "scen_II_221209/"),
                             model_wd = "L/", data_wd = "data/", rds = FALSE)

dat_L <- do.call(rbind, Map(cbind, it = seq_along(fit_L), 
                            lapply(fit_L, function(x) {
  x$predictions$mu_long[, c("id", "marker", "obstime", "Mean")]
})))

dat_L <- rbind(dat_L, fit_L[[1]]$simulations$mu_long %>%
                 mutate(Mean = mu, mu = NULL, it = 0)) %>%
  mutate(src = factor(ifelse(it == 0, "TRUE", "EST")))

ggplot(dat_L %>% arrange(desc(it), src),
       aes(x = obstime, y = Mean, group = it, color = src, 
                  alpha = src)) + 
  facet_grid(id ~ marker) +
  geom_line() +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme_bw() +
  labs(y = expression(mu(t)~","~hat(mu)(t)), alpha = NULL, color = NULL) +
  ggtitle("Different Fits based on Different Seeds and FPCs")

dat_L_obs <- do.call(rbind, Map(cbind, it = seq_along(fit_L), 
                                lapply(fit_L, function(x) {
  x$predictions$mu[, c("id", "marker", "obstime", "Mean")]
})))
dat_Lobs <- rbind(dat_L_obs, fit_L[[1]]$simulations$mu %>%
                    mutate(Mean = mu, mu = NULL, it = 0)) %>%
  mutate(src = factor(ifelse(it == 0, "TRUE", "EST")))

ggplot(dat_Lobs %>% arrange(desc(it), src),
       aes(x = obstime, y = Mean, group = it, color = src, 
           alpha = src)) + 
  facet_grid(id ~ marker) +
  geom_line() +
  scale_alpha_manual(values = c(0.1, 1)) +
  theme_bw() +
  labs(y = expression(mu(t)~","~hat(mu)(t)), alpha = NULL, color = NULL) +
  ggtitle("Different Fits based on Different Seeds and FPCs")


# Different Fits with Different Seeds For TRUE FPCs -------------------


fit_E <- sim_predict_longobs(obs = ids,
                             list.files(paste0(server_wd, "scen_II_221209/E/")), 
                             wd = paste0(server_wd, "scen_II_221209/"),
                             model_wd = "E/", data_wd = "data/", rds = FALSE)

dat_E <- do.call(rbind, Map(cbind, it = seq_along(fit_E), 
                            lapply(fit_E, function(x) {
                              x$predictions$mu_long[, c("id", "marker",
                                                        "obstime", "Mean")]
                            })))

dat_E <- rbind(dat_E, fit_E[[1]]$simulations$mu_long %>%
                 mutate(Mean = mu, mu = NULL, it = 0)) %>%
  mutate(src = factor(ifelse(it == 0, "TRUE", "EST")))


ggplot(dat_E, aes(x = obstime, y = Mean, group = it, color = src, 
                  alpha = src)) + 
  geom_line() +
  facet_grid(id ~ marker) +
  scale_alpha_manual(values = c(0.6, 1)) +
  theme_bw() +
  labs(y = expression(mu(t)~","~hat(mu)(t)), alpha = NULL, color = NULL) +
  ggtitle("Different Fits based on Different Seeds")



# Look at one model -------------------------------------------------------

dat_Lobs %>% filter(marker == "m1", id == "113", Mean < 0) %>% head()
m111 <- sim_predict_longobs_i(obs = ids, paste0("b", 111, ".Rdata"), 
                              wd = paste0(server_wd, "scen_II_221209/"),
                              model_wd = "L/", data_wd = "data/", rds = FALSE)
d111 <- m111$predictions$mu_long %>%
  left_join(m111$simulations$mu_long, by = c("id", "marker", "obstime")) %>%
  left_join(m111$simulations$mu %>% select(id, marker, obstime, y),
            by = c("id", "marker", "obstime") )

ggplot(d111, aes(x = obstime)) +
  facet_grid(id ~ marker) +
  geom_line(aes(y = Mean), linetype = "dashed") +
  geom_line(aes(y = mu)) +
  geom_point(aes(y = y)) +
  theme_bw()

# Is the MFPCA the same as is calculated in simulation
load(paste0(server_wd, "scen_II_221209/L/b111.Rdata"))
mfpc <- attr(b_est, "FPCs")
mfpc111 <- readRDS(paste0(server_wd, "../../fpc_estimation/simulation/",
                          "230124_parallel_to_simscenII/A_uncens/m111.rds"))
all.equal(mfpc, mfpc111)
plot(mfpc111$functions)
ggplot(b_est$model.frame %>% 
         pivot_longer(6:11), aes(x = obstime, y = value, group = id)) +
  geom_point() +
  facet_grid(name~marker) +
  theme_bw()
# looks like the first FPC might be switched?

fpc_dat <- readRDS(paste0(server_wd, "../../fpc_estimation/simulation/",
                          "230124_parallel_to_simscenII/eval_dat.rds"))
fpc_dat %>% filter(it == 111, type != "obs", setting == "uncens")
# e.g. angle is not an extreme angle


# compare to other model 110
m110 <- sim_predict_longobs_i(obs = ids, paste0("b", 110, ".Rdata"), 
                              wd = paste0(server_wd, "scen_II_221209/"),
                              model_wd = "L/", data_wd = "data/", rds = FALSE)
d110 <- m110$predictions$mu_long %>%
  left_join(m110$simulations$mu_long, by = c("id", "marker", "obstime")) %>%
  left_join(m110$simulations$mu %>% select(id, marker, obstime, y),
            by = c("id", "marker", "obstime") )

ggplot(d110, aes(x = obstime)) +
  facet_grid(id ~ marker) +
  geom_line(aes(y = Mean), linetype = "dashed") +
  geom_line(aes(y = mu)) +
  geom_point(aes(y = y)) +
  theme_bw()
mfpc110 <- readRDS(paste0(server_wd, "../../fpc_estimation/simulation/",
                          "230124_parallel_to_simscenII/A_uncens/m110.rds"))
plot(mfpc110$functions)



# What happens with attach function? --------------------------------------

load(paste0(server_wd, "scen_II_221209/data/d111.Rdata"))
sim_dat <- d_rirs$data 
d_rirs_est <- JMbamlss:::attach_wfpc(mfpc, d_rirs$data, n = 6)
str(d_rirs_est)
# apparently, the FPC functions are only added on to the data set
# -> Maybe that is where the problem arises?
# => Rerun the simulation: F2 and L2



# New Simulation Scenarios ------------------------------------------------

F100 <- sim_predict_longobs_i(obs = ids, paste0("b", 100, ".rds"), 
                              wd = paste0(server_wd, "scen_II_230117/"),
                              model_wd = "F2/", data_wd = "data/", rds = TRUE)
f100_dat <- F100$predictions$mu_long %>%
  left_join(F100$simulations$mu_long, by = c("id", "marker", "obstime")) %>%
  left_join(F100$simulations$mu %>% select(id, marker, obstime, y),
            by = c("id", "marker", "obstime") )

ggplot(f100_dat, aes(x = obstime)) +
  facet_grid(id ~ marker) +
  geom_line(aes(y = Mean), linetype = "dashed") +
  geom_line(aes(y = mu)) +
  geom_point(aes(y = y)) +
  theme_bw()
