
# PBC Data Analysis with JMbayes2 -----------------------------------------


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
library(JMbayes2)
library(JMbamlss)
library(tidyverse)

# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
# Logarithmic scale for all three markers
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years + sqrt(.Machine$double.eps),
         obstime = year + sqrt(.Machine$double.eps), marker = factor(marker),
         logy = log(y)) %>%
  select(id, survtime, event, sex, drug, age, marker, obstime, y, logy) %>%
  arrange(marker, id, obstime) %>%
  na.omit() %>%
  as.data.frame()

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(dplyr::n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis)

# Use logy for logarithmic transform
# Same dataset has to be used for all models and lme cannot handle missing
# values
p_long_jmb <- p_long %>%
  pivot_wider(id_cols = c(id, obstime, survtime, event, sex, age),
              values_from = logy, names_from = marker) %>%
  na.omit()

p_long_id <- p_long_jmb %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()


# Cox Model
CoxFit <- coxph(Surv(survtime, event) ~ age + sex, data = p_long_id)



# Random Intercept - Slope ------------------------------------------------

# Univariate longitudinal models
fm1 <- lme(serBilir ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
           random = ~ obstime | id)
fm2 <- lme(serChol ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
           random = ~ obstime | id)
fm3 <- lme(SGOT ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
           random = ~ obstime | id)


# the joint model that links all sub-models
jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "obstime",
               n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
saveRDS(jointFit, file = "inst/objects/pbc_jmb_rirs.Rds")


# Random Intercept - ns(slope) --------------------------------------------

# Univariate longitudinal models
fm1ns <- lme(serBilir ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
             random = ~ ns(obstime, df = 2) | id)
fm2ns <- lme(serChol ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
             random = ~ ns(obstime, df = 2) | id)
fm3ns <- lme(SGOT ~ ns(obstime, df = 3) + sex, data = p_long_jmb, 
             random = ~ ns(obstime, df = 2) | id)
# All models do not work for random effects with df = 3


# the joint model that links all sub-models
jointFit <- jm(CoxFit, list(fm1ns, fm2ns, fm3ns), time_var = "obstime",
               n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
saveRDS(jointFit, file = "inst/objects/pbc_jmb_ns.Rds")