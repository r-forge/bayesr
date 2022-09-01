library(JMbayes2)
library(tidyverse)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/C/sim110.Rdata"))
dat <- out$simdat$data %>% 
  pivot_wider(names_from = marker, values_from = y)# %>%
  #arrange(id, obstime)

coxm <- coxph(Surv(survtime, event) ~ x3, data = dat)
longm1 <- lme(m1 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))
longm2 <- lme(m2 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))
longm3 <- lme(m3 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))
longm4 <- lme(m4 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))
longm5 <- lme(m5 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))
longm6 <- lme(m6 ~ obstime + x3, random = ~ ns(obstime, 2) | id,
              data = dat, na.action = na.exclude,
              control = lmeControl(opt = "optim"))

mjm <- jm(coxm, list(longm1, longm2, longm3, longm4, longm5, longm6),
          time_var = "obstime")

library(JMbayes)
coxm <- coxph(Surv(survtime, event) ~ x3, data = dat, model = TRUE)
mixed <- mvglmer(list(m1 ~ obstime + x3 + (ns(obstime, 2) | id),
                      m2 ~ obstime + x3 + (ns(obstime, 2) | id),
                      m3 ~ obstime + x3 + (ns(obstime, 2) | id),
                      m4 ~ obstime + x3 + (ns(obstime, 2) | id),
                      m5 ~ obstime + x3 + (ns(obstime, 2) | id),
                      m6 ~ obstime + x3 + (ns(obstime, 2) | id)), 
                 data = dat,
                 families = rep(list(gaussian), 6))
mjm <- mvJointModelBayes(mvglmerObject = mixed, survObject = coxm, 
                         timeVar = "obstime")

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum bilirubin
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)

# a linear mixed model for the prothrombin time
fm2 <- lme(serChol ~ year * sex, data = pbc2, random = ~ year | id,
           na.action = na.omit)

# the joint model that links all sub-models
jointFit <- jm(CoxFit, list(fm1, fm2), time_var = "year",
               n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
