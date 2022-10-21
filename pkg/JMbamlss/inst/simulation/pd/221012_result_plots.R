library(tidyverse)

wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation/scen_I_130922")
load(paste0(wd, "/results.Rdata"))

r <- do.call(rbind, c(r_est, r_jmb, r_tru))


# Survival ----------------------------------------------------------------

ggplot(data = r %>% filter(predictor == "lambga", type == "Bias"),
       aes(y = value, x = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Lambda + Gamma")

ggplot(data = r %>% filter(predictor == "lambga", type == "MSE"),
       aes(y = value, x = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Lambda + Gamma")

ggplot(data = r %>% filter(predictor == "lambga", type == "Coverage" ),
       aes(y = value, x = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Lambda + Gamma")



# Alpha -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "alpha", type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  ggtitle("Bias: Alpha")

ggplot(data = r %>% filter(predictor == "alpha", type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 1.5)) +
  ggtitle("MSE: Alpha")

ggplot(data = r %>% filter(predictor == "alpha", type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Alpha")



# Mu ----------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "mu", type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Bias: Mu")

ggplot(data = r %>% filter(predictor == "mu", type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.075)) +
  ggtitle("MSE: Mu")

ggplot(data = r %>% filter(predictor == "mu", type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Mu")


# Bias by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Bias",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  ggtitle("Bias: Mu (by t)")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Bias",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(-0.05, 0.05)) +
  ggtitle("Bias: Mu (by t)")


# MSE by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(0, 0.1)) +
  ggtitle("MSE: Mu (by t)")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "MSE",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(0, 0.1)) +
  ggtitle("MSE: Mu (by t)")


# Coverage by t
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Coverage",
                           marker %in% c("m1", "m2", "m3"),
                           t %in% seq(0, 1, by = 0.1)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Coverage: Mu (by t)")
ggplot(data = r %>% filter(predictor == "mu_long",
                           type == "Coverage",
                           marker %in% c("m4", "m5", "m6"),
                           t %in% seq(0, 1, by = 0.05)),
       aes(y = value, x = t, col = model)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  ggtitle("Coverage: Mu (by t)")


# Sigma -------------------------------------------------------------------


ggplot(data = r %>% filter(predictor == "sigma", type == "Bias"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Bias: Sigma")

ggplot(data = r %>% filter(predictor == "sigma",type == "MSE"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("MSE: Sigma")

ggplot(data = r %>% filter(predictor == "sigma", type == "Coverage"),
       aes(y = value, x = marker, col = model)) +
  geom_boxplot() +
  ggtitle("Coverage: Sigma")




# -------------------------------------------------------------------------


# Other models ------------------------------------------------------------

library(tidyverse)

wd1 <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
             "share=volkmana.hub/JMbamlss/simulation/scen_I_051022")
load(paste0(wd1, "/results.Rdata"))
load()

r <- do.call(rbind, c(r_975, r_full, r_tru, r_jmb, r_est))


