library(tidyverse)

wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation")
result_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
                    "share=volkmana.hub/JMbamlss/simulation/",
                    "result_plots/scen_I/")
load(paste0(wd, "/results.Rdata"))

r <- do.call(rbind, c(r_est_95, r_jmb, r_tru_1, r_tru_95, r_tru_975)) %>%
  mutate(Model = factor(model, levels = c("jmb", "tru_1", "tru_975", "tru_95",
                                          "est_95"),
                        labels = c("JMB", "Tru_1", "Tru_975", "Tru_95",
                                   "Est_95")))

# Export as PDF (6 x 8 in)


# Survival ----------------------------------------------------------------

ggplot(data = r %>% filter(predictor == "lambga", type == "Bias"),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Lambda + Gamma") +
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(result_wd, "lambga_bias.png"), device = png)

ggplot(data = r %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Lambda + Gamma") +
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(result_wd, "lambga_mse.png"), device = png)

ggplot(data = r %>% filter(predictor == "lambga", type == "Coverage" ),
       aes(y = value, x = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Lambda + Gamma")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(result_wd, "lambga_coverage.png"), device = png)



# Alpha -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "alpha", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Alpha")+
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(result_wd, "alpha_bias.png"), device = png)

ggplot(data = r %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Alpha")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(result_wd, "alpha_mse.png"), device = png)

ggplot(data = r %>% filter(predictor == "alpha", type == "Coverage"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Alpha")
ggsave(paste0(result_wd, "alpha_coverage.png"), device = png)

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
ggsave(paste0(result_wd, "mu_bias.png"), device = png)

ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.075)) +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(result_wd, "mu_mse.png"), device = png)

ggplot(data = r %>% filter(predictor == "mu", type == "MSE",
                           model %in% c("jmb", "tru_1")),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(result_wd, "mu_mse_jmb_tru1.png"), device = png)

ggplot(data = r %>% filter(predictor == "mu", type == "Coverage"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Coverage: Mu")+
  ylab("Emp. Coverage") +
  theme_bw()
ggsave(paste0(result_wd, "mu_coverage.png"), device = png)


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
ggsave(paste0(result_wd, "mut_bias1.png"), device = png)
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
ggsave(paste0(result_wd, "mut_bias2.png"), device = png)


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
ggsave(paste0(result_wd, "mut_coverage1.png"), device = png)
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
ggsave(paste0(result_wd, "mut_coverage2.png"), device = png)


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
ggsave(paste0(result_wd, "mut_coverage1.png"), device = png)
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
ggsave(paste0(result_wd, "mut_coverage2.png"), device = png)


# Sigma -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "sigma", type == "Bias"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("Bias: Sigma")+
  ylab("Emp. Bias") +
  theme_bw()
ggsave(paste0(result_wd, "sigma_bias.png"), device = png)

ggplot(data = r %>% filter(predictor == "sigma",type == "MSE"),
       aes(y = value, x = marker, col = Model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")+
  ylab("Emp. MSE") +
  theme_bw()
ggsave(paste0(result_wd, "sigma_mse.png"), device = png)

ggplot(data = r %>% filter(predictor == "sigma", type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Sigma")
ggsave(paste0(result_wd, "sigma_coverage.png"), device = png)

r %>% 
  filter(predictor == "sigma", type == "Coverage") %>%
  group_by(Model, marker) %>%
  summarize(Mean = mean(value)) %>%
  pivot_wider(names_from = Model, values_from = Mean) %>%
  xtable::xtable(digits = 3) %>%
  xtable::print.xtable(include.rownames = FALSE)
