library(bamlss)
library(tidyverse)
library(manipulate)
# If slider doesn't work
# manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))



# Compare Longitudinal Fits -----------------------------------------------


load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/res_975.Rdata"))

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/data",
            "/d150.Rdata"))

dat <- rbind(data.frame(id = d_rirs$data$id,
                        obstime = d_rirs$data$obstime,
                        marker = d_rirs$data$marker,
                        mu = d_rirs$data$mu,
                        y = d_rirs$data$y,
                        y_hat = res_975[[1]]$mu$Mean,
                        type = "obs"),
             data.frame(id = d_rirs$data_full$id,
                        obstime = d_rirs$data_full$obstime,
                        marker = d_rirs$data_full$marker,
                        mu = d_rirs$data_full$mu,
                        y = NA,
                        y_hat = res_975[[1]]$mu_long$Mean,
                        type = "long"))


# Longitudinal prediction between mu and mu_t is identical
manipulate(ggplot(dat %>% filter(type == "obs", id == subj), aes(x = obstime)) +
             geom_point(data = dat %>% filter(type == "long", id == subj), 
                        aes(y = y_hat), col = "red") +
             geom_line(aes(y = mu)) +
             geom_line(aes(y = y_hat), col = "blue") +
             geom_point(aes(y = y), shape = 8) +
             facet_grid(marker~.),
           subj = slider(1, 150))



load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/res_full.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/res_est1.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/res_tru.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/res_jmb.Rdata"))

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/data",
            "/d102.Rdata"))

d_obs <- data.frame(id = d_rirs$data$id,
                    obstime = d_rirs$data$obstime,
                    marker = d_rirs$data$marker,
                    obs = d_rirs$data$y,
                    y = c(d_rirs$data$mu, res_jmb[[3]]$mu$Mean,
                          res_full[[3]]$mu$Mean, res_tru[[3]]$mu$Mean,
                          res_est[[3]]$mu$Mean),
                    type = rep(c("mu", "jmb", "b1", "b95", "best"), 
                               each = nrow(d_rirs$data)))

# Jmb and b1 are close to true longitudinal trajectory
manipulate(ggplot(d_obs %>% filter( id == subj), aes(x = obstime, y = y,
                                                     col = type)) +
             geom_line() +
             geom_point(aes(y = obs)) +
             facet_grid(marker~.),
           subj = slider(1, 150))

d_obs <- data.frame(id = d_rirs$data$id,
                    obstime = d_rirs$data$obstime,
                    marker = d_rirs$data$marker,
                    obs = d_rirs$data$mu,
                    y = c(res_jmb[[3]]$mu$Mean,
                          res_full[[3]]$mu$Mean, res_tru[[3]]$mu$Mean,
                          res_est[[3]]$mu$Mean),
                    type = rep(c("jmb", "b1", "b95", "best"), 
                               each = nrow(d_rirs$data))) %>%
  mutate(dif = obs - y)


ggplot(d_obs %>% filter(marker %in% c("m1", "m2", "m3"), 
                        obstime %in% seq(0.1, 1, by = 0.1)), 
       aes(x = as.factor(obstime), y = dif, col = type)) +
  geom_boxplot() +
  facet_grid(marker~., scales = "free")
ggplot(d_obs %>% filter(marker %in% c("m4", "m5", "m6"), 
                        obstime %in% seq(0.1, 1, by = 0.1)), 
       aes(x = as.factor(obstime), y = dif, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")

ggplot(d_obs %>% filter(marker %in% c("m1", "m2", "m3"), 
                        obstime %in% seq(0.1, 1, by = 0.1),
                        type %in% c("b1", "jmb")), 
       aes(x = as.factor(obstime), y = dif, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")
ggplot(d_obs %>% filter(marker %in% c("m4", "m5", "m6"), 
                        obstime %in% seq(0.1, 1, by = 0.1),
                        type %in% c("b1", "jmb")),  
       aes(x = as.factor(obstime), y = dif, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")



load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/results.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/results.Rdata"))

d_obs %>% 
  filter(type == "b1") %>%
  group_by(marker) %>% 
  summarise(Mean = mean(dif)) %>%
  ungroup()
r_full[[3]] %>% 
  filter(type == "Bias", predictor == "mu")


d_long <- data.frame(id = d_rirs$data_full$id,
                    obstime = d_rirs$data_full$obstime,
                    marker = d_rirs$data_full$marker,
                    obs = d_rirs$data_full$mu,
                    y = c(res_jmb[[3]]$mu_long$Mean,
                          res_full[[3]]$mu_long$Mean, res_tru[[3]]$mu_long$Mean,
                          res_est[[3]]$mu_long$Mean),
                    type = rep(c("jmb", "b1", "b95", "best"), 
                               each = nrow(d_rirs$data_full))) %>%
  mutate(dif = obs - y)


loo <- r_full[[3]] %>% 
  filter(type == "Bias", predictor == "mu_long") %>%
  select(value)
man <- d_long %>% 
  filter(type == "b1") %>%
  group_by(obstime, marker) %>%
  summarise(Mean = mean(dif)) %>%
  ungroup() %>%
  select(Mean) %>%
  mutate(Mean = -Mean)
all.equal(loo, man)

   data.frame(id = d_rirs$data_full$id,
              obstime = d_rirs$data_full$obstime,
              marker = d_rirs$data_full$marker,
              mu = d_rirs$data_full$mu,
              y = NA,
              y_hat = res_975[[1]]$mu_long$Mean,
              type = "long"))


# Reconstruct the calculation ---------------------------------------------

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/res_full.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/sim_dat_full.Rdata"))
source("R/sim_helperfun.R")

r_b1 <- sim_results(res_full, sim_dat_full, name = "b1")

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022",
            "/results.Rdata"))
all.equal(r_b1, r_full)

##############
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/data",
            "/d102.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/res_jmb.Rdata"))
d_long <- data.frame(id = d_rirs$data_full$id,
                     obstime = d_rirs$data_full$obstime,
                     marker = d_rirs$data_full$marker,
                     obs = d_rirs$data_full$mu,
                     y = c(res_jmb[[3]]$mu_long$Mean,
                           res_full[[3]]$mu_long$Mean),
                     type = rep(c("jmb", "b1"), 
                                each = nrow(d_rirs$data_full))) %>%
  mutate(squ_err = (obs - y)^2)
ggplot(d_long %>% filter(marker %in% c("m1", "m2", "m3"), 
                        obstime %in% seq(0.0, 1, by = 0.1)), 
       aes(x = as.factor(obstime), y = squ_err, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")
ggplot(d_long %>% filter(marker %in% c("m4", "m5", "m6"), 
                        obstime %in% seq(0, 1, by = 0.1)),  
       aes(x = as.factor(obstime), y = squ_err, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")

d_obs <- data.frame(id = d_rirs$data$id,
                    obstime = d_rirs$data$obstime,
                    marker = d_rirs$data$marker,
                    obs = d_rirs$data$mu,
                    y = c(res_jmb[[3]]$mu$Mean,
                          res_full[[3]]$mu$Mean),
                    type = rep(c("jmb", "b1"), 
                               each = nrow(d_rirs$data))) %>%
  mutate(squ_err = (obs - y)^2)
ggplot(d_obs %>% filter(marker %in% c("m1", "m2", "m3")), 
       aes(x = as.factor(obstime), y = squ_err, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")
ggplot(d_obs %>% filter(marker %in% c("m4", "m5", "m6")),  
       aes(x = as.factor(obstime), y = squ_err, col = type)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(marker~., scales = "free")
d_obs %>% group_by(obstime, marker) %>% summarize(n = n()) %>% ungroup()

# Is it only better when looking at the actual observed times
aha <- left_join(d_long, d_obs, by = c("id", "marker", "obstime", "type"))
a <- aha %>% 
  mutate(diff = squ_err.x - squ_err.y)
sum(a$diff, na.rm = TRUE)
a %>% filter(diff != 0) %>% head(n= 20)
# DIFFERENT PREDICTIONS!