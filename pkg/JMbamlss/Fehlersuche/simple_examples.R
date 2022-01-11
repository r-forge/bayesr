
# Different Simple Data Examples ------------------------------------------

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
#debugonce(simMultiJM)


# Univariate JMs ----------------------------------------------------------

set.seed(1808)

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
var(dat1$data_full$mu) # variation in longitudinal trajectories 
var(dat1$data_full$s1)/120 # comes from PCRE
ggplot(dat1$data_hypo, aes(x = obstime, y = mu, group = id)) + 
  geom_line() # constant trajectories


# Constant high hazard -> everyone dies in the beginning
dat1_1 <- simMultiJM(nmark = 1, M = 5,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
ggplot(dat1_1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_point() # all die at 0
var(dat1_1$data_full$mu)
var(dat1_1$data_full$s1) # no longer the same
var(dat1_1$data_full$s5) # eigenvalues are seq(1, 0.2, by = 0.2)
ggplot(dat1$data_hypo, aes(x = obstime, y = mu, group = id)) + 
  geom_line() # constant trajectories


# Constant high hazard -> everyone dies in the beginning but different time 
# scale
dat1_2 <- simMultiJM(nmark = 1, M = 1,
                     times = seq(0, 1, length.out = 121),
                     lambda = function(t, x) 300,
                     gamma = function(x) 0,
                     alpha = list(function(t, x) 0*t),
                     mu = list(function(t, x) 1.25),
                     sigma = function(t, x) 0.001 + 0*t,
                     full = TRUE)
table(dat1_2$data$event) # no censoring
ggplot(dat1_2$data, aes(x = obstime, y = mu, group = id)) + 
  geom_point() # all die at 0
var(dat1_2$data_full$mu)
var(dat1_2$data_full$s1) # variation comes from PCRE


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
table(dat2_1$data$survtime < 120) # (almost) no one is censored
ggplot(dat2_1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))
ggplot(dat2_1$data, aes(x = obstime, y = mu, group = id)) + 
  geom_point() # missingness shares of longitudinal trajectories
summary(as.integer(table(dat2_1$data$id))) # 0.25*121 = 30.25



# Fit for simple model ----------------------------------------------------

set.seed(1808)

# Constant hazard
dat3 <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                   lambda = function(t, x) -5,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
ggplot(dat3$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

# debug(preproc_MFPCA)
# MFPCA on the data
mfpca1 <- preproc_MFPCA(dat3$data)
mfpcaT <- create_true_MFPCA(M = 1, nmarker = 1)

f_pcre <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime +
    s(id, wfpc, bs = "pcre", xt = list("mfpc" = mfpcaT)),
  sigma ~ 1,
  alpha ~ 1
)
b_pcre <- bamlss(f_pcre, family = mjm_bamlss, data = dat3$data, 
                 timevar = "obstime", sampler = FALSE, maxit = 200)

f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
b_re <- bamlss(f_re, family = mjm_bamlss, data = dat3$data, timevar = "obstime",
               sampler = FALSE, maxit = 200)

f_old <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
b_old <- bamlss(f_old, data = dat3$data, family = "jm", timevar = "obstime",
                idvar = "id", sampler = FALSE, maxit = 200)

# Obvious differences
matrix(c(b_pcre$parameters$lambda$s[[1]], b_re$parameters$lambda$s[[1]],
         b_old$parameters$lambda$s[[1]]), nrow = 3, byrow = TRUE)
plot(b_re, model = "lambda", ask = FALSE)
plot(b_old, model = "lambda", ask = FALSE)


##  Only survival part
# Fix the longitudinal part
b_nolong <- bamlss(f_pcre, family = mjm_bamlss, data = dat3$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200, 
                   opt_long = FALSE)
matrix(c(b_pcre$parameters$lambda$s[[1]], b_nolong$parameters$lambda$s[[1]]),
       nrow = 2, byrow = TRUE)

b_re_nolong <- bamlss(f_re, family = mjm_bamlss, data = dat3$data, 
                      timevar = "obstime", sampler = FALSE, maxit = 200,
                      opt_long = FALSE)

b_old_fix <- bamlss(f_old, data = dat3$data, family = "jmFEHLERSUCHE",
                    timevar = "obstime",
                idvar = "id", sampler = FALSE, maxit = 200, fix.alpha = TRUE,
                fix.mu = TRUE, fix.dalpha = TRUE, fix.sigma = TRUE)

matrix(c(b_old$parameters$lambda$s[[1]], b_old_fix$parameters$lambda$s[[1]]),
       nrow = 2, byrow = TRUE)
plot(b_re_nolong, model = "lambda", ask = FALSE)
plot(b_old, model = "lambda", ask = FALSE)
plot(b_old_fix, model = "lambda", ask = FALSE)



## Only survival part with bad integral
# Fix the longitudinal part and reduce subdivisions
b_old_fix_bad <- bamlss(f_old, data = dat3$data, family = "jmFEHLERSUCHE",
                        timevar = "obstime", idvar = "id", sampler = FALSE,
                        maxit = 200, fix.alpha = TRUE, fix.mu = TRUE, 
                        fix.dalpha = TRUE, fix.sigma = TRUE, subdivisions = 7)
b_old_fix_bad$parameters$lambda$s[[1]]


## Only survival part with fixed smoothing parameter in bamlss
# Fix the longitudinal part and fix the smoothing parameter
b_old_fix_smo <- bamlss(f_old, data = dat3$data, family = "jmFEHLERSUCHE",
                        timevar = "obstime", idvar = "id", sampler = FALSE,
                        maxit = 200, fix.alpha = TRUE, fix.mu = TRUE, 
                        fix.dalpha = TRUE, fix.sigma = TRUE, do.optim2 = FALSE)
b_old_fix_smo$parameters$lambda
plot(b_old_fix_smo, model = "lambda", ask = FALSE)



# Change the intercept in the formulas ------------------------------------

f_new_wint <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
b_new_wint <- bamlss(f_new_wint, family = mjm_bamlss, data = dat3$data,
                     timevar = "obstime", sampler = FALSE, maxit = 200, 
                     opt_long = FALSE)
matrix(c(b_new_wint$parameters$lambda$s[[1]], 
         b_re_nolong$parameters$lambda$s[[1]]), nrow = 2, byrow = TRUE)

f_old_wint <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
b_old_wint <-  bamlss(f_old_wint, data = dat3$data, family = "jmFEHLERSUCHE",
                      timevar = "obstime", idvar = "id", sampler = FALSE, 
                      maxit = 200, fix.alpha = TRUE, fix.mu = TRUE, 
                      fix.dalpha = TRUE, fix.sigma = TRUE)
matrix(c(b_old_wint$parameters$lambda$s[[1]], 
         b_old_fix$parameters$lambda$s[[1]]), nrow = 2, byrow = TRUE)

b_old_wint_nosmoo <- bamlss(f_old_wint, data = dat3$data, family = "jmFEHLERSUCHE",
                            timevar = "obstime", idvar = "id", sampler = FALSE,
                            maxit = 200, fix.alpha = TRUE, fix.mu = TRUE, 
                            fix.dalpha = TRUE, fix.sigma = TRUE,
                            do.optim2 = FALSE)
matrix(c(b_old_wint_nosmoo$parameters$lambda$s[[1]], 
         b_old_fix_smo$parameters$lambda$s[[1]]), nrow = 2, byrow = TRUE)


# Fit for example with linear baseline ------------------------------------

set.seed(1808)

# Constant hazard
dat4 <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                   lambda = function(t, x) -5 + 0.02*t,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
ggplot(dat4$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

## Linear specification of baseline hazard doesn't work in bamlss implementation
# flin_re <- list(
#   Surv2(survtime, event, obs = y) ~ survtime,
#   gamma ~ 1,
#   mu ~ 1 + obstime + s(id, bs = "re"),
#   sigma ~ 1,
#   alpha ~ 1
# )
# blin_re <- bamlss(flin_re, family = mjm_bamlss, data = dat4$data,
#                   timevar = "obstime", sampler = FALSE, maxit = 200, 
#                   opt_long = FALSE)
# 
# flin_old <- list(
#   Surv2(survtime, event, obs = y) ~ -1 + survtime,
#   gamma ~ 1,
#   mu ~ 1 + obstime + s(id, bs = "re"),
#   sigma ~ 1,
#   alpha ~ 1,
#   dalpha ~ -1
# )
# blin_old <-  bamlss(flin_old, data = dat3$data, family = "jmFEHLERSUCHE",
#                       timevar = "obstime", idvar = "id", sampler = FALSE, 
#                       maxit = 200, fix.alpha = TRUE, fix.mu = TRUE, 
#                       fix.dalpha = TRUE, fix.sigma = TRUE)

b_lin_re <- bamlss(f_re, family = mjm_bamlss, data = dat4$data,
                   timevar = "obstime", sampler = FALSE, maxit = 200,
                   opt_long = FALSE)
b_lin_old <- bamlss(f_old, data = dat4$data, family = "jmFEHLERSUCHE",
                    timevar = "obstime", idvar = "id", sampler = FALSE,
                    maxit = 200, fix.alpha = TRUE, fix.mu = TRUE,
                    fix.dalpha = TRUE, fix.sigma = TRUE)
b_lin_old_fixsm <- bamlss(f_old, data = dat4$data, family = "jmFEHLERSUCHE",
                          timevar = "obstime", idvar = "id", sampler = FALSE,
                          maxit = 200, fix.alpha = TRUE, fix.mu = TRUE,
                          fix.dalpha = TRUE, fix.sigma = TRUE,
                          do.optim2 = FALSE)



# Try Cox Model in bamlss -------------------------------------------------


set.seed(1808)

# Constant hazard
dat5 <- simMultiJM(nsub = 3000, nmark = 1, M = 1,
                   lambda = function(t, x) -5 + 0.02*t,
                   gamma = function(x) 0,
                   alpha = list(function(t, x) 0*t),
                   mu = list(function(t, x) 1.25),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
ggplot(dat5$data, aes(x = obstime, y = mu, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

f_cox <- list(
  Surv(survtime, event) ~ s(survtime),
  gamma ~ 1
)
dat_cox <- dat5$data %>%
  distinct(id, event, survtime)
b_cox <- bamlss(f_cox, data = dat_cox, family = "cox", maxit = 100,
                sampler = FALSE)
plot(b_cox)
# Plot shows that there is a linear baseline hazard with rough slope of 
# 2.2/120 = 0.0183
b_cox$parameters$gamma



# Multivariate JMs --------------------------------------------------------



dat2 <- simMultiJM(nmark = 2, M = 3,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = rep(list(function(t, x) 0*t), 2),
                   mu = rep(list(function(t, x) 1.25), 2),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)

dat3 <- simMultiJM(nmark = 3, M = 1,
                   lambda = function(t, x) 300,
                   gamma = function(x) 0,
                   alpha = rep(list(function(t, x) 0*t), 3),
                   mu = rep(list(function(t, x) 1.25), 3),
                   sigma = function(t, x) 0.001 + 0*t,
                   full = TRUE)
