library(bamlss)
library(tidyverse)
library(manipulate)
# If slider doesn't work
# manipulate(plot(1:5, cex=size), size = slider(0.5,10,step=0.5))
source("R/preprocessing.R")
source("R/eval_mfun.R")



load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/data",
            "/d102.Rdata"))
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/bamlss_tru/b102.Rdata"))

# Predict based on model frame
pred_modframe <- predict(b_est, 
                         newdata = b_est$model.frame,
                         model = "mu")


# Predict based on simulated data
mfpca_est <- attr(b_est, "FPCs")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
newdat_orig <- na.omit(attach_wfpc(mfpca_est, d_rirs$data, n = nfpc))
pred_origdat <- predict(b_est,
                        newdata = newdat_orig,
                        model = "mu")
all.equal(pred_modframe, pred_origdat)
newdat_orig <- newdat_orig %>%
  mutate(yhat = pred_origdat)


# Predict based on full data
newdat_full <- na.omit(attach_wfpc(mfpca_est, d_rirs$data_full, n = nfpc))
pred_fulldat <- predict(b_est,
                        newdata = newdat_full,
                        model = "mu")
newdat_full <- newdat_full %>% 
  mutate(ytil = pred_fulldat)
comp <- right_join(newdat_full, newdat_orig, by = c("id", "marker", "obstime"))
all.equal(comp$mu.x, comp$mu.y)
all.equal(comp$fpc.1.x, comp$fpc.1.y)
all.equal(comp$ytil, comp$yhat)

newdat_orig <- newdat_orig %>%
  mutate(mse_orig = (mu - yhat)^2)
mean(newdat_orig$mse_orig)

newdat_full <- newdat_full %>%
  mutate(mse_full = (mu - ytil)^2)
mean(newdat_full$mse_full)



load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922",
            "/results.Rdata"))
r_man_full <- newdat_full %>% group_by(obstime, marker) %>%
  summarise(MEAN = mean(mse_full)) %>%
  ungroup()
r_cal <- r_tru[[3]] %>%
  filter(type == "MSE", model == "tru", predictor == "mu_long")
all.equal(r_man_full$MEAN, r_cal$value)

r_man_orig <- newdat_orig %>% group_by(obstime, marker) %>%
  summarise(MEAN = mean(mse_orig)) %>%
  ungroup()

ggplot(rbind(r_man_full, r_man_orig) %>% 
         mutate(dat = factor(rep(c("full", "orig"), each = 606))), 
       aes(x = obstime, y = MEAN, color = dat)) +
  geom_boxplot() +
  facet_wrap(marker ~ ., scale = "free")


