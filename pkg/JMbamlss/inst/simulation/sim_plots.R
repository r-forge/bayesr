
# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


# Always
library(survival)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/MJM_opt.R")
source("R/opt_updating_cpp.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing_cpp.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/compile.R")
compile_alex(location)


ids <-  c(23, 45, 67)
set.seed(1808)
n_ids <- sample(1:150, 4)

# 6 Dimensions but more wiggly estimate
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/B/sim110.Rdata"))
dat <- out$simdat$data_full %>%
  mutate(fit_t = predict(out$b_tru, newdata = out$simdat$data_full, 
                         model = "mu"),
         fit_e = predict(out$b_est, newdata = out$simdat$data_full, 
                         model = "mu"))
ggplot(dat %>% filter(id %in% ids), aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu), linetype = "dotted") +
  geom_line(aes(y = fit_t), linetype = "dashed") +
  geom_line(aes(y = fit_e)) +
  facet_wrap(~ marker, scales = "free") +
  geom_point(data = out$simdat$data %>% filter(id %in% ids), aes(y = y))
ggplot(dat %>% filter(id %in% n_ids), aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu), linetype = "dotted") +
  geom_line(aes(y = fit_t), linetype = "dashed") +
  geom_line(aes(y = fit_e)) +
  facet_wrap(~ marker, scales = "free") +
  geom_point(data = out$simdat$data %>% filter(id %in% n_ids), aes(y = y))


# 6 Dimensions
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/C/sim110.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/C/jmbayes110.Rdata"))
dat <- out$simdat$data_full %>%
  mutate(#fit_t = predict(out$b_tru, newdata = out$simdat$data_full, 
          #               model = "mu"),
         fit_e = predict(out$b_est, newdata = out$simdat$data_full, 
                         model = "mu"))
ggplot(dat %>% filter(id %in% ids), aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu), linetype = "dotted") +
  #geom_line(aes(y = fit_t), linetype = "dashed") +
  geom_line(aes(y = fit_e)) +
  facet_wrap(~ marker, scales = "free") +
  geom_point(data = out$simdat$data %>% filter(id %in% ids), aes(y = y))
ggplot(dat %>% filter(id %in% n_ids), aes(x = obstime, colour = id)) +
  geom_line(aes(y = mu), linetype = "dotted") +
  #geom_line(aes(y = fit_t), linetype = "dashed") +
  geom_line(aes(y = fit_e)) +
  facet_wrap(~ marker, scales = "free") +
  geom_point(data = out$simdat$data %>% filter(id %in% n_ids), aes(y = y))
aha <- predict(mjm, newdata = out$simdat$data_full)
