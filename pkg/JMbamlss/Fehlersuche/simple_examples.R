
# Different Simple Data Examples ------------------------------------------

# Data Generation
source("R/simMultiJM.R")

# Helperfunction PCRE
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")

library(tidyverse)
debugonce(simMultiJM)



# Univariate JMs ----------------------------------------------------------


# Constant high hazard -> everyone dies in the beginning
dat1 <- simMultiJM(nmark = 1, M = 1,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
table(dat1$data$event) # no censoring
ggplot(dat1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_point() # all die at 0
sd(dat1$data_full$mu)
sd(dat1$data_full$s1) # variation in longitudinal trajectores comes from PCRE
ggplot(dat1$data_hypo, aes(x = obstime, y = mu, group = id)) + 
  geom_line() # constant trajectories


# Constant low hazard -> everyone survives
dat2 <- simMultiJM(nmark = 1, M = 1,
                   lambda = function(t, x) -300,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
table(dat2$data$event) # no one dies
table(dat2$data$survtime < 120) # half are censored
ggplot(dat2$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))


# Constant low hazard with (almost) no censoring
dat2_1 <- simMultiJM(nmark = 1, M = 1, maxfac = 500,
                     lambda = function(t, x) -300,
                     gamma = function(x) 0,
                     alpha = list(function(t, x) 0*t),
                     mu = list(function(t, x) 1.25),
                     sigma = function(t, x) 0.001 + 0*t,
                     full = TRUE)
table(dat2_1$data$event) # no one dies
table(dat2_1$data$survtime < 120) # no one is censored
ggplot(dat2_1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))
ggplot(dat2_1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_point() # missingness shares of longitudinal trajectories
summary(as.integer(table(dat2_1$data$id))) # 0.25*121 = 30.25






# Multivariate JMs --------------------------------------------------------



dat2 <- simMultiJM(nmark = 2, M = 3,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = rep(list(function(t, x) 0*t), 2),
                   mu = rep(list(function(t, x) 1.25), 2),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)

dat3 <- simMultiJM(nmark = 2, M = 1,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = rep(list(function(t, x) 0*t), 2),
                   mu = rep(list(function(t, x) 1.25), 2),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
