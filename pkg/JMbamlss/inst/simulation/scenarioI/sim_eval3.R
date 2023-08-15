

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


# Evaluate bamlss TRUE FPCs -----------------------------------------------

mnames_btru <- list.files(path = file.path(server_wd, "scen_I_230719", 
                                           "bamlss_tru"))
preds_btru <- JMbamlss:::sim_bamlss_predict(mnames_btru, server_wd, 
                                            "/scen_I_230719/bamlss_tru/",
                                            "/scen_I_230719/data/", rds = TRUE,
                                            old = TRUE)
saveRDS(preds_btru, 
        file = file.path(server_wd, "scen_I_230719", "preds_btru.rds"))

it_list <- JMbamlss:::sim_results(lapply(preds_btru, "[[", "predictions"),
                                  lapply(preds_btru, "[[", "simulations"),
                                  name = "TRU")
eval_btru <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
saveRDS(eval_btru,
        file = file.path(server_wd, "scen_I_230719", "eval_btru.rds"))
rm(preds_btru, eval_btru, it_list)



# Evaluate bamlss EST 1 FPCs ----------------------------------------------

mnames_best1 <- list.files(path = file.path(server_wd, "scen_I_230719", 
                                            "bamlss_est1"))
preds_best1 <- JMbamlss:::sim_bamlss_predict(mnames_best1, server_wd, 
                                             "/scen_I_230719/bamlss_est1/",
                                             "/scen_I_230719/data/", rds = TRUE,
                                             old = TRUE)
saveRDS(preds_best1, 
        file = file.path(server_wd, "scen_I_230719", "preds_best1.rds"))

it_list <- JMbamlss:::sim_results(lapply(preds_best1, "[[", "predictions"),
                                  lapply(preds_best1, "[[", "simulations"),
                                  name = "EST1")
eval_best1 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                 it_list))
saveRDS(eval_best1,
        file = file.path(server_wd, "scen_I_230719", "eval_best1.rds"))
rm(preds_best1, eval_best1, it_list)




# Evaluate bamlss EST 95 FPCs ---------------------------------------------

mnames_best95 <- list.files(path = file.path(server_wd, "scen_I_230719", 
                                            "bamlss_est95"))
preds_best95 <- JMbamlss:::sim_bamlss_predict(mnames_best95, server_wd, 
                                             "/scen_I_230719/bamlss_est95/",
                                             "/scen_I_230719/data/", rds = TRUE,
                                             old = TRUE)
saveRDS(preds_best95, 
        file = file.path(server_wd, "scen_I_230719", "preds_best95.rds"))

it_list <- JMbamlss:::sim_results(lapply(preds_best95, "[[", "predictions"),
                                  lapply(preds_best95, "[[", "simulations"),
                                  name = "EST95")
eval_best95 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                 it_list))
saveRDS(eval_best95,
        file = file.path(server_wd, "scen_I_230719", "eval_best95.rds"))
rm(preds_best95, eval_best95, it_list)



# Evaluate jmbayes  -------------------------------------------------------

mnames_jmb <- list.files(path = file.path(server_wd, "scen_I_230719", "jmb"))
preds_jmb <- JMbamlss:::sim_jmb_predict(mnames_jmb, server_wd, 
                                        "/scen_I_230719/jmb/",
                                        "/scen_I_230719/data/", rds = TRUE)
saveRDS(preds_jmb, 
        file = file.path(server_wd, "scen_I_230719", "preds_jmb.rds"))

it_list <- JMbamlss:::sim_results(lapply(preds_jmb, "[[", "predictions"),
                                  lapply(preds_jmb, "[[", "simulations"),
                                  name = "JMB")
eval_jmb <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
saveRDS(eval_jmb,
        file = file.path(server_wd, "scen_I_230719", "eval_jmb.rds"))
rm(preds_jmb, eval_jmb, it_list)



# Table for Paper ---------------------------------------------------------

e_btru <- readRDS(file.path(server_wd, "scen_I_230719", "eval_btru.rds"))
e_best1 <- readRDS(file.path(server_wd, "scen_I_230719", "eval_best1.rds"))
e_best95 <- readRDS(file.path(server_wd, "scen_I_230719", "eval_best95.rds"))
e_jmb <- readRDS(file.path(server_wd, "scen_I_230719", "eval_jmb.rds"))


r <- rbind(e_btru, e_best1, e_best95, e_jmb) %>%
  mutate(Model = factor(model, levels = c("TRU", "EST1", "EST95", "JMB"),
                        labels = c("TRUE", "EST", "TRUNC", "JMB")))

bias <- r %>%
  filter(type == "Bias", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  mutate(mean = ifelse(predictor == "mu", mean*1000, mean)) %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  rowid_to_column() %>%
  mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
         marker = ifelse(predictor == "sigma", "all", marker)) %>%
  group_by(rowid, predictor, marker) %>%
  summarise(across(1:4, mean), .groups = "drop") %>%
  slice(7, 1:6, 8:14)
  
mse <- r %>%
  filter(type == "MSE", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  mutate(mean = ifelse(predictor == "mu", mean*1000, mean)) %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  rowid_to_column() %>%
  mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
         marker = ifelse(predictor == "sigma", "all", marker)) %>%
  group_by(rowid, predictor, marker) %>%
  summarise(across(1:4, mean), .groups = "drop") %>%
  slice(7, 1:6, 8:14)

coverage <- r %>%
  filter(type == "Coverage", 
         !(predictor %in% c("mu_long", "lambga", "lambga_long"))) %>%
  group_by(Model, predictor, marker) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  pivot_wider(id_cols = c(predictor, marker), names_from = Model, 
              values_from = mean) %>%
  rowid_to_column() %>%
  mutate(rowid = ifelse(predictor == "sigma", 16, rowid),
         marker = ifelse(predictor == "sigma", "all", marker)) %>%
  group_by(rowid, predictor, marker) %>%
  summarise(across(1:4, mean), .groups = "drop") %>%
  slice(7, 1:6, 8:14)

eval <- left_join(bias, mse, by = c("rowid", "predictor", "marker"), 
                  suffix = c("bias", "mse")) %>%
  left_join(coverage, by = c("rowid", "predictor", "marker")) %>%
  mutate(rowid = NULL,
         marker = NULL,
         predictor = c("$\\lambda + \\gamma$", paste0("$\\alpha_{", 1:6, "}$"),
                       paste0("$\\mu_{", 1:6, "}$"), "$\\sigma$"))

xtable::print.xtable(xtable::xtable(eval, digits = 3), include.rownames = FALSE,
                     sanitize.text.function = identity)



# Plot longitudinal fits --------------------------------------------------

d_rirs <- readRDS(file.path(server_wd, "scen_I_230719", "data", "d100.rds"))
set.seed(1444)
ids <- c(sample(which(d_rirs$data_short$survtime[1:150] > 0.2 & 
                        d_rirs$data_short$survtime[1:150] < 0.35), size = 1),
         sample(which(d_rirs$data_short$survtime[1:150] > 0.35 & 
                        d_rirs$data_short$survtime[1:150] < 0.7), size = 1),
         sample(which(d_rirs$data_short$survtime[1:150] > 0.7), size = 1))

p_btru <- readRDS(file.path(server_wd, "scen_I_230719", 
                            "preds_btru.rds"))$b100.rds$predictions
p_best1 <- readRDS(file.path(server_wd, "scen_I_230719", 
                             "preds_best1.rds"))$b100.rds$predictions
p_best95 <- readRDS(file.path(server_wd, "scen_I_230719", 
                              "preds_best95.rds"))$b100.rds$predictions
p_jmb <- readRDS(file.path(server_wd, "scen_I_230719", 
                           "preds_jmb.rds"))$jmb100.rds$predictions

# Compare the fitted values to the observed values
mean_dat <- d_rirs$data_full %>%
  select(id, marker, obstime, mu) %>%
  left_join(d_rirs$data %>% select(id, marker, obstime, y), 
            by = c("id", "marker", "obstime")) %>%
  cbind(p_btru$mu_long %>% mutate("TRUE" = Mean) %>%
          select("TRUE")) %>%
  cbind(p_best1$mu_long %>% mutate(EST = Mean) %>% 
          select(EST)) %>%
  cbind(p_best95$mu_long %>% mutate(TRUNC = Mean) %>%
          select(TRUNC)) %>%
  cbind(p_jmb$mu_long %>% mutate(JMB = Mean) %>%
          select(JMB)) %>%
  pivot_longer(c(4, 6:ncol(.))) %>%
  mutate(name = factor(name, levels = c("mu", "TRUE", "EST", "TRUNC", "JMB")))

ggplot(mean_dat %>% filter(id %in% ids) %>%
         mutate(marker = factor(marker, labels = paste("Outcome", 1:6)),
                id = factor(id, labels = c("Subject 28", "Subject 70", 
                                           "Subject 130"))),
       aes(x = obstime, y = value, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = c("black", scales::hue_pal()(4))) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL,
       x = "Time")
# save 6x12
