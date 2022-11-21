library(tidyverse)

wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation")
plot_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berli",
                  "n.de,share=volkmana.hub/JMbamlss/simulation/",
                    "result_plots/scen_I_est/")
load(paste0(wd, "/results.Rdata"))
load(paste0(wd, "/results_est.Rdata"))

r <- do.call(rbind, c(r_est_95, r_est_975, r_est_1, r_jmb)) %>%
  mutate(Model = factor(model, levels = c("jmb", "est_1", "est_975", "est_95"),
                        labels = c("JMB", "Est_1", "Est_975", "Est_95")))

# Export as PDF (6 x 8 in)


# Survival ----------------------------------------------------------------

ggplot(data = r %>% filter(predictor == "lambga", type == "Bias"),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Lambda + Gamma") +
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(plot_wd, "lambga_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Lambda + Gamma") +
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "lambga_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "lambga", type == "Coverage" ),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Lambda + Gamma")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(plot_wd, "lambga_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")



# Alpha -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "alpha", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Alpha")+
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(plot_wd, "alpha_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Alpha")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "alpha_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "alpha", type == "Coverage"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Alpha")
ggsave(paste0(plot_wd, "alpha_coverage.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

r %>% 
  filter(predictor == "alpha", type == "Coverage") %>%
  group_by(Model, marker) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = Model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)


# Mu ----------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "mu", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Bias: Mu")+
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "mu", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Bias: Mu (Outliers removed)")+
  ylab("Emp. Bias") +
  scale_y_continuous(limits = c(-0.01, 0.01)) + 
  theme_bw()
ggsave(paste0(plot_wd, "mu_bias_outliers.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_mse.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.075)) +
  ggtitle("MSE: Mu (Outliers removed)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_mse_outliers.pdf"), device = "pdf", width = 8, 
       height = 6, units = "in")

ggplot(data = r %>% filter(predictor == "mu", type == "Coverage"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Mu")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(plot_wd, "mu_coverage.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")


# Bias by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Bias",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  ggtitle("Bias: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_bias1.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Bias",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  ggtitle("Bias: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_bias2.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")


# MSE by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(0, 0.1)) +
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
  scale_y_continuous(limits = c(0, 0.1)) +
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_mse2.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")


# Coverage by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Coverage",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Coverage: Mu (by t)")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_coverage1.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Coverage",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.05)),
       aes(y = value, x = t, col = Model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Coverage: Mu (by t)")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(plot_wd, "mut_coverage2.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")



# Sigma -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "sigma", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Bias: Sigma")+
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(plot_wd, "sigma_bias.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "sigma",type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(plot_wd, "sigma_mse.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

ggplot(data = r %>% filter(predictor == "sigma", type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Sigma")
ggsave(paste0(plot_wd, "sigma_coverage.pdf"), device = "pdf", width = 8, height = 6,
       units = "in")

r %>% 
  filter(predictor == "sigma", type == "Coverage") %>%
  group_by(Model, marker) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = Model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)
