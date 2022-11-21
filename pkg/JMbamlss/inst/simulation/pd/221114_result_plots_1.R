location <- "workstation"

if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if (location == "server_linux") {
    paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/simulation/")
  } else NULL
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
  plot_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                    "berlin.de,share=volkmana.hub/JMbamlss/simulation/",
                    "result_plots/scen_I_tru1/")
}

library(survival)
library(MFPCA)
library(tidyverse)
library(devtools)
load_all()

# Prediction Part ---------------------------------------------------------

m_tru_1 <- list.files(path = paste0(results_wd, "/scen_I_051022/bamlss_tru_1"))

pred_tru_1 <- sim_bamlss_predict(m_tru_1, results_wd, 
                                 m_setting = "/scen_I_051022/", 
                                 d_setting = "/scen_I_130922/", 
                                 folder = "bamlss_tru_1/")

# Simulated Data Part -----------------------------------------------------

d_tru_1 <- paste0("d", substr(m_tru_1, 2, 10))

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


# Comparing Predictions to Data -------------------------------------------

r_tru_1_new <- sim_results(pred_tru_1, sdat_tru_1, name = "tru_1_new")
save(r_tru_1_new, file = paste0(results_wd, "/results_tru_1.Rdata"))


# Comparing Results for different Codes -----------------------------------

load(paste0(results_wd, "/results.Rdata"))
load(paste0(results_wd, "/results_est.Rdata"))
r <- do.call(rbind, c(r_jmb, r_tru_1, r_tru_1_new)) %>%
  mutate(Model = factor(model, levels = c("jmb", "tru_1", "tru_1_new"),
                        labels = c("JMB", "Tru_1 Old", "Tru_1 New")))
r_comb <- do.call(rbind, c(r_tru_1, r_tru_1_new)) %>%
  mutate(Model = factor(model, levels = c("tru_1", "tru_1_new"),
                        labels = c("Tru_1", "Tru_1")))
r <- rbind(r, r_comb) %>%
  mutate(Model = factor(Model, levels = c("JMB", "Tru_1", "Tru_1 Old", 
                                          "Tru_1 New")))


# Lambga
ggplot(data = r %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Lambda + Gamma") +
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "lambga_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")


# Alpha
ggplot(data = r %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Alpha")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "alpha_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

# Mu
ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
      aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_mse.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

# Mu by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse1.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse2.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

# Sigma
ggplot(data = r %>% filter(predictor == "sigma",type == "MSE"),
       aes(y = log(value), x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")+
  ylab("Log(Emp. MSE)") +
  theme_bw()
ggsave(paste0(plot_wd, "sigma_mse.pdf"), device = "pdf", width = 8,
       height = 6, units = "in")

