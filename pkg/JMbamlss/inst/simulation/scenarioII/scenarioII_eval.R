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



# Evaluate the simulation -------------------------------------------------


m_bamlss <- JMbamlss:::sim_jmbamlss_eval(
  wd = paste0(server_wd, "scen_II_221209/"), model_wd = "A/", data_wd = "data/",
  name = "JMbamlss")
save(m_bamlss, file = paste0(server_wd, "scen_II_221209/res_A.Rdata"))

# m_jmbayes <- JMbamlss:::sim_jmbayes_eval(
#   wd = paste0(server_wd, "scen_II_221209/"), model_wd = "A_jmb/", 
#   data_wd = "data/", name = "JMbayes")
# save(m_jmbayes, file = paste0(server_wd, "scen_II_221209/res_A_jmb.Rdata"))

m_jmbayes <- JMbamlss:::sim_jmbayes_eval(
  wd = paste0(server_wd, "scen_II_221209/"), model_wd = "B_jmb/", 
  data_wd = "data/", name = "JMbayes")
save(m_jmbayes, file = paste0(server_wd, "scen_II_221209/res_B_jmb.Rdata"))


# Plot the results --------------------------------------------------------

load(paste0(server_wd, "scen_II_221209/res_A.Rdata"))
load(paste0(server_wd, "scen_II_221209/res_A_jmb.Rdata"))

r <- rbind(m_bamlss, m_jmbayes)

# Lambga
ggplot(data = r %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = model)) +
  geom_boxplot() +
  ggtitle("MSE: Lambda + Gamma") +
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "lambga_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

# Alpha
ggplot(data = r %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Alpha")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "alpha_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

# Mu
ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_mse.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

# Mu by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           t %in% seq(0, 1, by = 0.05)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

# Sigma
ggplot(data = r %>% filter(predictor == "sigma", type == "MSE"),
       aes(y = log(value), x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")+
  ylab("Log(Emp. MSE)") +
  theme_bw()
ggsave(paste0(plot_wd, "sigma_mse.pdf"), device = "pdf", width = 8,
       height = 6, units = "in")
