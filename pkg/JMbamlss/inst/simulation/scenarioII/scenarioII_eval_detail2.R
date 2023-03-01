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

f100 <- JMbamlss:::sim_bamlss_predict_i("b100.rds", 
                                        paste0(server_wd, "scen_II_230117/"),
                                        "F2/", "data/", rds = TRUE)
j100 <- JMbamlss:::sim_jmb_predict_i("jmb100.rds", 
                                     paste0(server_wd, "scen_II_230117/"),
                                     "JMB/", "data/", rds = TRUE)
l100 <- JMbamlss:::sim_bamlss_predict_i("b100.rds", 
                                        paste0(server_wd, "scen_II_230117/"),
                                        "L2/", "data/", rds = TRUE)
e100 <- JMbamlss:::sim_bamlss_predict_i("b100.rds", 
                                        paste0(server_wd, "scen_II_230117/"),
                                        "E/", "data/", rds = TRUE)
a100 <- JMbamlss:::sim_bamlss_predict_i("b100.rds", 
                                        paste0(server_wd, "scen_II_230117/"),
                                        "A2/", "data/", rds = TRUE)
d_sim <- readRDS(paste0(server_wd, "scen_II_230117/data/d100.rds"))
m_a100 <- readRDS(paste0(server_wd, "scen_II_230117/A2/b100.rds"))

# Compare the fitted values to the observed values
mean_dat <- d_sim$data_full %>%
  select(id, marker, obstime, mu) %>%
  left_join(d_sim$data %>% select(id, marker, obstime, y), 
            by = c("id", "marker", "obstime")) %>%
  cbind(f100$predictions$mu_long %>% mutate(CenFPC1 = Mean) %>%
          select(CenFPC1)) %>%
  cbind(j100$predictions$mu_long %>% mutate(JMB = Mean) %>% 
          select(JMB)) %>%
  cbind(l100$predictions$mu_long %>% mutate(UncenFPC1 = Mean) %>%
          select(UncenFPC1)) %>%
  cbind(e100$predictions$mu_long %>% mutate(TruFPC1 = Mean) %>%
          select(TruFPC1)) %>%
  cbind(a100$predictions$mu_long %>% mutate(CenFPC95 = Mean) %>%
          select(CenFPC95)) %>%
  pivot_longer(c(4, 6:ncol(.))) %>%
  mutate(name = factor(name, levels = c("mu", "TruFPC1", "CenFPC1",
                                        "UncenFPC1", "CenFPC95", "JMB")))


# Random observations
set.seed(1312)
ids <- sample(1:300, size = 5)

# Plot all fits together
ggplot(mean_dat %>% filter(id %in% ids),
       aes(x = obstime, y = value, linetype = name, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "purple", "red", "blue", "orange",
                                "green")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL) +
  ggtitle("Simulation Scenario II: Results for Mu",
          "Random Longitudinal Trajectories and Fits on Simulation Run 1")


# Plot different FPC sources
ggplot(mean_dat %>% filter(id %in% ids, name %in% c("mu", "TruFPC1", "CenFPC1", 
                                                    "UncenFPC1")),
       aes(x = obstime, y = value, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "purple", "red", "blue")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL) +
  ggtitle("Simulation Scenario II: Results for Mu",
          "Random Longitudinal Trajectories and Fits on Simulation Run 1")

# Plot different basis functions
ggplot(mean_dat %>% filter(id %in% ids, name %in% c("mu", "CenFPC1", "CenFPC95",
                                                    "JMB")),
       aes(x = obstime, y = value, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "red", "orange", "green")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL) +
  ggtitle("Simulation Scenario II: Results for Mu",
          "Random Longitudinal Trajectories and Fits on Simulation Run 1")

# For Finland
ggplot(mean_dat %>% filter(id %in% c(56, 147), 
                           name %in% c("mu", "CenFPC1", "JMB")) %>%
         mutate(name = factor(name, labels = c("Truth", "Proposed",
                                               "Other")),
                marker = factor(marker, labels = c("Exposure 1", "Exposure 2")),
                id = factor(id, labels = c("Subject 56", "Subject 147"))),
       aes(x = obstime, y = value, color = name)) +
  facet_grid(id ~ marker, scales = "free_y")+
  theme_bw() +
  geom_point(aes(y = y), color = "grey") +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "blue", "red")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL,
       x = "Time") +
  ggtitle("Examples of Model Fits to Simulated Data")
# save as (3x9inch pdf)



# Evaluate the simulation settings ----------------------------------------

scenII_A2 <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "A2/", data_wd = "data/", name = "A", rds = TRUE)
saveRDS(scenII_A2, file = paste0(server_wd, "scen_II_230117/res_A2.rds"))

scenII_F2 <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "F2/", data_wd = "data/", name = "F", rds = TRUE)
saveRDS(scenII_F2, file = paste0(server_wd, "scen_II_230117/res_F2.rds"))

scenII_L2 <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "L2/", data_wd = "data/", name = "L", rds = TRUE)
saveRDS(scenII_L2, file = paste0(server_wd, "scen_II_230117/res_L2.rds"))

scenII_E <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "E/", data_wd = "data/", name = "E", rds = TRUE)
saveRDS(scenII_E, file = paste0(server_wd, "scen_II_230117/res_E.rds"))

scenII_JMB <- JMbamlss:::sim_jmbayes_eval(
  wd = paste0(server_wd, "../../JMbamlss/simulation/scen_II_230117/"),
  model_wd = "JMB/", data_wd = "data/", name = "JMB", rds = TRUE)
saveRDS(scenII_JMB, file = paste0(server_wd, "scen_II_230117/res_JMB.rds"))

scenII_A2 <- readRDS(paste0(server_wd, "scen_II_230117/res_A2.rds"))
scenII_F2 <- readRDS(paste0(server_wd, "scen_II_230117/res_F2.rds"))
scenII_L2 <- readRDS(paste0(server_wd, "scen_II_230117/res_L2.rds"))
scenII_E <- readRDS(paste0(server_wd, "scen_II_230117/res_E.rds"))
scenII_JMB <- readRDS(paste0(server_wd, "scen_II_230117/res_JMB.rds"))

eval_dat <- rbind(scenII_A2, scenII_F2, scenII_L2, scenII_E, scenII_JMB) %>%
  mutate(model = factor(model, levels = c("E", "L", "F", "A", "JMB"),
                        labels = c("True", "Uncens", "Cens", "Cens95", "JMB")))

# Lambda + Gamma predictor ---
ggplot(data = eval_dat %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = model, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("MSE: Lambda + Gamma")
ggsave(paste0(plot_wd, "lambga_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% filter(predictor == "lambga", type == "Bias"),
       aes(y = value, x = model, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. Bias") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("Bias: Lambda + Gamma")
ggsave(paste0(plot_wd, "lambga_bias.pdf"), device = "pdf", width = 8,
       height = 6, units = "in")

ggplot(data = eval_dat %>% filter(predictor == "lambga", type == "Coverage" ),
       aes(y = value, x = model, color = model)) +
  geom_boxplot() +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  theme(legend.position = "None") +
  ggtitle("Coverage: Lambda + Gamma")
ggsave(paste0(plot_wd, "lambga_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

# Alpha predictor ---
ggplot(data = eval_dat %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +  
  theme_bw() + 
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("MSE: Alpha")
ggsave(paste0(plot_wd, "alpha_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% filter(predictor == "alpha", type == "Bias"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Emp. Bias") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ggtitle("Bias: Alpha")
ggsave(paste0(plot_wd, "alpha_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

eval_dat %>% 
  filter(predictor == "alpha", type == "Coverage") %>%
  group_by(model, marker) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)

# Mu predictor ---
ggplot(data = eval_dat %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave(paste0(plot_wd, "mu_mse.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

ggplot(data = eval_dat %>% filter(predictor == "mu", type == "Bias"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +
  ggtitle("Bias: Mu")+
  ylab("Emp. Bias") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave(paste0(plot_wd, "mu_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% filter(predictor == "mu", type == "Coverage"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Mu")+
  ylab("Emp. Coverage") +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dotted")
ggsave(paste0(plot_wd, "mu_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

# Mu by t predictor ---
ggplot(data = eval_dat %>% 
         filter(predictor == "mu_long",
                type == "MSE",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")
# cutoff
ggplot(data = eval_dat %>% 
         filter(predictor == "mu_long",
                type == "MSE",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("MSE: Mu (by t) - Cutoff at 0.2")+
  ylab("Emp. MSE") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 0.2)) +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse_cutoff.pdf"), device = "pdf", width = 8,
       height = 6, units = "in")

ggplot(data = eval_dat %>% 
         filter(predictor == "mu_long",
                type == "Bias",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Bias: Mu (by t)")+
  ylab("Emp. Bias") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% 
         filter(predictor == "mu_long",
                type == "Coverage",
                t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Coverage: Mu (by t)")+
  ylab("Emp. Coverage") +
  geom_hline(yintercept = 0.95, linetype = "dotted") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")


# Sigma ---
ggplot(data = eval_dat %>% filter(predictor == "sigma", type == "MSE"),
       aes(y = value, x = marker, color = model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")+
  ylab("Emp. MSE") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave(paste0(plot_wd, "sigma_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = eval_dat %>% filter(predictor == "sigma", type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Bias: Sigma")+
  ylab("Emp. Bias") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dotted")
ggsave(paste0(plot_wd, "sigma_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

eval_dat %>% 
  filter(predictor == "sigma", type == "Coverage") %>%
  group_by(model, marker) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)

