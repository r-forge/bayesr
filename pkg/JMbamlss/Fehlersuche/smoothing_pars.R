
# Analyze influence of smoothing parameters -------------------------------


# Data Generation and Preprocessing
source("R/simMultiJM.R")
source("R/preprocessing.R")

# Helperfunction PCRE
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")

# Alternative bamlss Code (more options)
source("Fehlersuche/JM.R")

library(tidyverse)
library(refund)
library(survival)
library(bamlss)



# Simple model with linear baseline hazard --------------------------------

set.seed(1808)

# Linear hazard, no association, constant marker values
dat0 <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                   lambda = function(t, x) -5 + 0.02*t,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
ggplot(dat0$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))


f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)

b_re <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, timevar = "obstime",
               sampler = FALSE, maxit = 200, opt_long = FALSE)
b_re$parameters$lambda$s$`s(survtime)`[10]

f_sp <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
b_sp <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, timevar = "obstime",
               sampler = FALSE, maxit = 200, opt_long = FALSE, 
               tau = list("lambda" = list("s(survtime)" = 0.005)))
b_sp2 <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, timevar = "obstime",
                sampler = FALSE, maxit = 200, opt_long = FALSE, 
                tau = list("lambda" = list("s(survtime)" = 5000)))

plot(b_re, ask = FALSE, model = "lambda")
plot(b_sp, ask = FALSE, model = "lambda")
plot(b_sp2, ask = FALSE, model = "lambda")
# If we approximate with linear function, then range:
# b_re:  (-2, 1) -- 3/120 = 0.025
# b_sp:  (-1, 1.5) -- 2.5/120 = 0.021
# b_sp2: (-3, 3) -- 6/120 = 0.05
# First two are close to the true slope of 0.02
c(b_re$parameters$gamma[[1]], b_sp$parameters$gamma[[1]],
  b_sp2$parameters$gamma[[1]])


b_jm <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, timevar = "obstime",
               sampler = FALSE, maxit = 200)
b_jm_sp <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, 
                  timevar = "obstime", sampler = FALSE, maxit = 200,
                  tau = list("lambda" = list("s(survtime)" = 0.005)))
b_jm_sp2 <- bamlss(f_re, family = mjm_bamlss, data = dat0$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200,
                   tau = list("lambda" = list("s(survtime)" = 5000)))
plot(b_jm, ask = FALSE, model = "lambda")
plot(b_jm_sp, ask = FALSE, model = "lambda")
plot(b_jm_sp2, ask = FALSE, model = "lambda")
# b_re:  (-2.5, 2) -- 4.5/120 = 0.038
# b_sp:  (-1.5, 2.5) -- 4/120 = 0.033
# b_sp2: (-4, 4) -- 8/120 = 0.07
c(b_jm$parameters$gamma[[1]], b_jm_sp$parameters$gamma[[1]],
  b_jm_sp2$parameters$gamma[[1]])
# Model cannot be identified
c(b_jm$parameters$alpha[[1]], b_jm_sp$parameters$alpha[[1]],
  b_jm_sp2$parameters$alpha[[1]])


set.seed(1808)

# Linear hazard, no association, constant marker values
dat1 <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                   lambda = function(t, x) -5 + 0.02*t,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 100 + 0*t,
                   full = TRUE)
ggplot(dat1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

