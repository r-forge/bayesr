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

m100 <- JMbamlss:::sim_bamlss_predict_i("b100.Rdata", 
                                        paste0(server_wd, "scen_II_221209/"),
                                        "A/", "data/")
j100 <- JMbamlss:::sim_jmb_predict_i("jmb100.Rdata", 
                                     paste0(server_wd, "scen_II_221209/"),
                                     "A_jmb/", "data/")
j100b <- JMbamlss:::sim_jmb_predict_i("jmb100.Rdata", 
                                     paste0(server_wd, "scen_II_221209/"),
                                     "B_jmb/", "data/")
load(paste0(server_wd, "scen_II_221209/data/d100.Rdata"))
load(paste0(server_wd, "scen_II_221209/A/b100.Rdata"))

set.seed(1312)
ids <- sample(1:300, size = 3)

# Compare the fitted values to the observed values
mean_dat <- d_rirs$data %>%
  select(id, marker, obstime, mu, y) %>%
  cbind(m100$predictions$mu %>% select(Mean)) %>%
  cbind(j100$predictions$mu %>% mutate(j_Mean = Mean) %>% select(j_Mean))
ggplot(mean_dat %>% filter(id %in% ids),
       aes(x = obstime, color = id)) +
  facet_wrap(~marker) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(y = y)) +
  geom_line(aes(y = mu)) +
  geom_point(aes(y = Mean), shape = 3) +
  geom_point(aes(y = j_Mean), shape = 1)

# Compare the fitted values to the underlying true values
mean_t_dat <- d_rirs$data_full %>%
  select(id, marker, obstime, mu) %>%
  cbind(m100$predictions$mu_long %>% select(Mean)) %>%
  cbind(j100$predictions$mu_long %>% mutate(j_Mean = Mean) %>% select(j_Mean))
ggplot(mean_t_dat %>% filter(id %in% ids),
       aes(x = obstime, color = id)) +
  facet_wrap(~marker) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_line(aes(y = mu)) +
  geom_point(aes(y = Mean), shape = 3) +
  geom_point(aes(y = j_Mean), shape = 1)

set.seed(1520)
ids <- sample(1:300, size = 5)
ggplot(mean_t_dat %>% filter(id %in% ids) %>%
         pivot_longer(4:6) %>%
         mutate(name = factor(name, 
                              levels = c("mu", "Mean", "j_Mean"),
                              labels = c("Data", "JMbamlss", 
                                         "JMbayes"))),
       aes(x = obstime, y = value, linetype = name, color = name)) +
  facet_grid(id~marker, scales = "free_y")+
  theme_bw() +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "red", "blue")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL) +
  ggtitle("Simulation Scenario II: Results for Mu",
          "Random Longitudinal Trajectories and Fits")

library(manipulate)
# manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))
manipulate({ggplot(mean_t_dat %>% filter(id == i) %>%
                     pivot_longer(4:6) %>%
                     mutate(name = factor(name, 
                                          levels = c("mu", "Mean", "j_Mean"),
                                          labels = c("Data", "JMbamlss", 
                                                     "JMbayes"))),
       aes(x = obstime, y = value, linetype = name, color = name)) +
  facet_wrap(~marker)+
  theme_bw() +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "red", "blue")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL)},
  i = slider(1, 300))



# Add up the effects for id 108 -------------------------------------------


# Predictions seem to work just fine
id108 <- d_rirs$data_full %>% filter(id == 108)
muhat108 <- predict(b_est, newdata = id108, model = "mu")
mufpc108 <- predict(b_est, newdata = id108, model = "mu", term = "id")
fpc1_108 <- predict(b_est, newdata = id108, model = "mu", term = "fpc.1")
fpc2_108 <- predict(b_est, newdata = id108, model = "mu", term = "fpc.2")
fpc3_108 <- predict(b_est, newdata = id108, model = "mu", term = "fpc.3")
fpc4_108 <- predict(b_est, newdata = id108, model = "mu", term = "fpc.4")

all.equal(fpc1_108 + fpc2_108 + fpc3_108 + fpc4_108,
          mufpc108)
id_108 <- id108 %>%
  cbind(muhat = muhat108, fpchat = mufpc108, fpc1hat = fpc1_108, 
        fpc2hat = fpc2_108, fpc3hat = fpc3_108, fpc4hat = fpc4_108) %>%
  mutate(cum0 = muhat - fpchat, cum1 = cum0 + fpc1hat, cum2 = cum1 + fpc2hat,
         cum3 = cum2 + fpc3hat, cum4 = cum3 + fpc4hat)%>%
  select(id, marker, obstime, mu, 35:45) %>%
  pivot_longer(5:15)
ggplot(id_108, aes(x = obstime, y = value)) +
  facet_grid(marker~name, scales = "free_y")+
  theme_bw() +
  geom_line() +
  scale_x_continuous(limits = c(0, 1))



# Evaluate both models on simulation 100 ----------------------------------

mse <- d_rirs$data %>%
  select(id, marker, obstime) %>%
  mutate(mse_jmbamlss = (m100$simulations$mu$mu - m100$predictions$mu$Mean)^2,
         mse_jmbayes = (j100$simulations$mu$mu - j100$predictions$mu$Mean)^2) %>%
  pivot_longer(4:5)
mse_long <- d_rirs$data_full %>%
  select(id, marker, obstime) %>%
  mutate(mse_jmbamlss = (m100$simulations$mu_long - 
                           m100$predictions$mu_long$Mean)^2,
         mse_jmbayes = (j100$simulations$mu_long - 
                          j100$predictions$mu_long$Mean)^2) %>%
  pivot_longer(4:5)


ggplot(data = mse,
       aes(y = value, x = marker, col = name)) +
  geom_boxplot() +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()

ggplot(data = mse_long,
       aes(y = value, x = marker, col = name)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("MSE: Mu")+
  ylab("Emp. MSE") +
  theme_bw()

ggplot(data = mse_long %>% filter(obstime %in% seq(0, 1, by = 0.05)),
       aes(y = value, x = factor(obstime), col = name)) +
  geom_boxplot() +
  facet_grid(marker ~., scales = "free") +
  scale_y_continuous(limits = c(0, 10))+
  ggtitle("MSE: Mu (by t)")+
  ylab("Emp. MSE") +
  theme_bw()


# Are the predictions even the same?
mse_joint <- mse %>%
  left_join(mse_long, by = c("id", "marker", "obstime", "name")) %>%
  mutate(diff = value.x - value.y)
head(mse_joint[which(mse_joint$diff != 0), ], n = 20)


# For JMBayes, the predictions for the data differ from the predictions for the
# simulated data
pred_short <- d_rirs$data %>% 
  select(id, marker, obstime, mu) %>%
  cbind(pred_short = j100$predictions$mu$Mean)
pred_long <- d_rirs$data_full %>%
  select(id, marker, obstime, mu) %>%
  cbind(pred_long = j100$predictions$mu_long$Mean)
pred <- pred_short %>%
  left_join(pred_long, by = c("id", "marker", "obstime", "mu")) %>%
  mutate(diff = pred_short - pred_long)
head(pred[which(pred$diff != 0), ])

# Create indicator which rows of the design matrices for the full are needed in
# the actual data

available <- d_rirs$data_full %>%
  select(id, marker, obstime) %>%
  left_join(d_rirs$data %>% select(id, marker, obstime, mu),
            by = c("id", "marker", "obstime"))
av <- !is.na(available$mu)

debug(JMbamlss:::sim_jmb_predict_i)
j100 <- JMbamlss:::sim_jmb_predict_i("jmb100.Rdata", 
                                     paste0(server_wd, "scen_II_221209/"),
                                     "A_jmb/", "data/")
all.equal(X[[1]], X_long[[1]][av[1:18487], ]) # only attributes differ
all.equal(X[[2]], X_long[[2]][av[18488:36974], ]) # only attributes differ
all.equal(Z[[1]], Z_long[[1]][av[1:18487], ])
all.equal(Z[[2]], Z_long[[2]][av[18488:36974], ]) # HERE THE RANDOM EFFECTS 
# DESIGN MATRICES DIFFER!
# Reason: Misspecification in the JMbayes fit: Original fit use bs() on m1 but
# ns() on m2
# Refit the model in B_jmb: Then no problemo

mse_b <- d_rirs$data %>%
  select(id, marker, obstime) %>%
  mutate(mse_jmbamlss = (m100$simulations$mu$mu - m100$predictions$mu$Mean)^2,
         mse_jmbayes = (j100$simulations$mu$mu - 
                          j100b$predictions$mu$Mean)^2) %>%
  pivot_longer(4:5)
mse_long_b <- d_rirs$data_full %>%
  select(id, marker, obstime) %>%
  mutate(mse_jmbamlss = (m100$simulations$mu_long - 
                           m100$predictions$mu_long$Mean)^2,
         mse_jmbayes = (j100$simulations$mu_long - 
                          j100b$predictions$mu_long$Mean)^2) %>%
  pivot_longer(4:5)

# Are the predictions even the same?
mse_joint_b <- mse_b %>%
  left_join(mse_long_b, by = c("id", "marker", "obstime", "name")) %>%
  mutate(diff = value.x - value.y)
head(mse_joint_b[which(mse_joint$diff != 0), ], n = 20)
# YES, NOW THEY ARE THE SAME!!



# Better Fit for Higher Cutoff? -------------------------------------------


m100b <- JMbamlss:::sim_bamlss_predict_i("b100.Rdata", 
                                         paste0(server_wd, "scen_II_221209/"),
                                         "B/", "data/")
m100c <- JMbamlss:::sim_bamlss_predict_i("b100.Rdata", 
                                         paste0(server_wd, "scen_II_221209/"),
                                         "C/", "data/")
m100d <- JMbamlss:::sim_bamlss_predict_i("b100.Rdata", 
                                         paste0(server_wd, "scen_II_221209/"),
                                         "D/", "data/")
load(paste0(server_wd, "scen_II_221209/A/b100.Rdata"))
b_est_A <- b_est
load(paste0(server_wd, "scen_II_221209/B/b100.Rdata"))
b_est_B <- b_est
load(paste0(server_wd, "scen_II_221209/C/b100.Rdata"))

set.seed(1520)
ids <- sample(1:300, size = 5)

mean_t_dat <- d_rirs$data_full %>%
  select(id, marker, obstime, mu) %>%
  cbind(m100$predictions$mu_long %>% mutate(Bam95 = Mean) %>% select(Bam95)) %>%
  cbind(m100b$predictions$mu_long %>% mutate(Bam975 = Mean) %>% select(Bam975)) %>%
  cbind(m100c$predictions$mu_long %>% mutate(Bam99 = Mean) %>% select(Bam99)) %>%
  cbind(m100d$predictions$mu_long %>% mutate(Bam1 = Mean) %>% select(Bam1)) %>%
  cbind(j100b$predictions$mu_long %>% mutate(JMB = Mean) %>% select(JMB)) %>%
  mutate(Mu = mu)

ggplot(mean_t_dat %>% filter(id %in% ids) %>%
         pivot_longer(5:10) %>%
         mutate(name = factor(name, 
                              levels = c("Mu", "Bam95", "Bam975", "Bam99", 
                                         "Bam1", "JMB"))),
       aes(x = obstime, y = value, linetype = name, color = name)) +
  facet_grid(id~marker, scales = "free_y")+
  theme_bw() +
  geom_line() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("black", "red", "orange", "purple", "blue", "green")) +
  labs(y = expression(mu(t)~","~hat(mu)(t)), linetype = NULL, color = NULL) +
  ggtitle("Simulation Scenario II: Results for Mu",
          "Random Longitudinal Trajectories and Fits")
