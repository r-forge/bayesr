
# Analyze influence of smoothing parameters -------------------------------
# UNIVARIATE

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

library(Matrix)
library(tidyverse)
library(refund)
library(survival)
library(bamlss)
library(funData)



# Linear BH, Small Sigma --------------------------------------------------

set.seed(1808)

# Linear hazard, no association, constant marker values
dat_sigmas <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                         lambda = function(t, x) -5 + 0.02*t,
                         gamma = function(x) 0,
                         alpha = list(function(t, x) 0*t),
                         mu = list(function(t, x) 1.25),
                         sigma = function(t, x) -50 + 0*t,
                         full = TRUE)
ggplot(dat_sigmas$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

# Create data for prediction
pred_dat <- dat_sigmas$data[c(1, 1), c("id", "survtime", "obstime")]
pred_dat$survtime <- pred_dat$obstime <- c(0, 120)

f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)

# Different models with different smoothing parameters
b_nolong <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200, 
                   opt_long = FALSE)
b_nolong$parameters$lambda$s$`s(survtime)`[10]

b_nolong_spl <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       opt_long = FALSE, 
                       tau = list("lambda" = list("s(survtime)" = 0.005)))
b_nolong_spu <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       opt_long = FALSE, 
                       tau = list("lambda" = list("s(survtime)" = 5000)))

# Compare models
plot(b_nolong, ask = FALSE, model = "lambda")
plot(b_nolong_spl, ask = FALSE, model = "lambda")
plot(b_nolong_spu, ask = FALSE, model = "lambda")
# If we approximate with linear function, then range:
# b_nolong:  (-2, 1) -- 3/120 = 0.025
# b_nolong_spl:  (-1, 1.5) -- 2.5/120 = 0.021
# b_nolong_spu: (-3, 3) -- 6/120 = 0.05
# First two are close to the true slope of 0.02

c(b_nolong$parameters$gamma[[1]], b_nolong_spl$parameters$gamma[[1]],
  b_nolong_spu$parameters$gamma[[1]])

# True hazard: [-5, -2.6]
# b_nolong_spl very close (but arbitrary smoothing parameter!)
predict(b_nolong, pred_dat)$lambda + predict(b_nolong, pred_dat)$gamma
predict(b_nolong_spl, pred_dat)$lambda + predict(b_nolong_spl, pred_dat)$gamma
predict(b_nolong_spu, pred_dat)$lambda + predict(b_nolong_spu, pred_dat)$gamma




# Comparison to bamlss ----------------------------------------------------


f_rebamlss <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
bamlss_nolong <- bamlss(f_rebamlss, family = "jmFEHLERSUCHE", 
                        data = dat_sigmas$data, timevar = "obstime", 
                        idvar = "id", sampler = FALSE, maxit = 200, 
                        fix.alpha = TRUE, fix.mu = TRUE, fix.sigma = TRUE)
(sp_bamlss <- bamlss_nolong$parameters$lambda$s$`s(survtime)`[10])

b_nolong_bamls <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                         timevar = "obstime", sampler = FALSE, maxit = 200,
                         opt_long = FALSE, 
                         tau = list("lambda" = list("s(survtime)" = sp_bamlss)))
plot(bamlss_nolong, ask = FALSE, model = "lambda")
plot(b_nolong_bamls, ask = FALSE, model = "lambda")
# bamlss_nolong: (-0.15, 0.25) -- 0.4/120 = 0.003
# b_nolong_bamls: (-0.01, 0.02) -- 0.03/120 = 0.00025
bamlss_nolong$parameters$gamma$p

# True hazard: [-5, -2.6]
predict(bamlss_nolong, pred_dat)$lambda + 
  predict(bamlss_nolong, pred_dat)$gamma
predict(b_nolong_bamls, pred_dat)$lambda + 
  predict(b_nolong_bamls, pred_dat)$gamma




# Include Longitudinal Data -----------------------------------------------


b_long <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                 timevar = "obstime", sampler = FALSE, maxit = 200,)
b_long$parameters$lambda$s$`s(survtime)`[10] # Same smoothing parameter

b_long_spl <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data,
                     timevar = "obstime", sampler = FALSE, maxit = 200,
                     tau = list("lambda" = list("s(survtime)" = 0.005)))
b_long_spu <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data,
                     timevar = "obstime", sampler = FALSE, maxit = 200,
                     tau = list("lambda" = list("s(survtime)" = 5000)))

# Estimated BH does not change
plot(b_long, ask = FALSE, model = "lambda")
plot(b_long_spl, ask = FALSE, model = "lambda")
plot(b_long_spu, ask = FALSE, model = "lambda")

# Identifiability problems with constant marker values and gamma intercept
c(b_long$parameters$gamma[[1]], b_long_spl$parameters$gamma[[1]],
  b_long_spu$parameters$gamma[[1]])
c(b_long$parameters$alpha[[1]], b_long_spl$parameters$alpha[[1]],
  b_long_spu$parameters$alpha[[1]])
b_long$parameters$gamma[[1]] + 
  b_long$parameters$alpha[[1]]*b_long$parameters$mu$p[[1]]

# True hazard: [-5, -2.6]
predict(b_long, pred_dat)$lambda + 
  predict(b_long, pred_dat)$gamma +
  predict(b_long, pred_dat)$alpha * predict(b_long, pred_dat)$mu
predict(b_long_spl, pred_dat)$lambda + 
  predict(b_long_spl, pred_dat)$gamma +
  predict(b_long_spl, pred_dat)$alpha * predict(b_long_spl, pred_dat)$mu
predict(b_long_spu, pred_dat)$lambda + 
  predict(b_long_spu, pred_dat)$gamma +
  predict(b_long_spu, pred_dat)$alpha * predict(b_long_spu, pred_dat)$mu



# Center Longitudinal Data ------------------------------------------------

dat_cen <- dat_sigmas
dat_cen$data$y <- dat_cen$data$y - mean(dat_cen$data$y)

b_long_cen <- bamlss(f_re, family = mjm_bamlss, data = dat_cen$data, 
                     timevar = "obstime", sampler = FALSE, maxit = 200)

# Estimated BH changes considerably
# (-8, 8) -- 16/120 = 0.13
plot(b_long_cen, ask = FALSE, model = "lambda")

# Alpha part is highly biased in own implementation
# But overall the hazard is close to the true hazard
b_long_cen$parameters$gamma$p
b_long_cen$parameters$alpha$p
b_long_cen$parameters$mu$p
b_long_cen$parameters$gamma[[1]] + 
  b_long_cen$parameters$alpha[[1]]*b_long_cen$parameters$mu$p[[1]]

# True hazard: [-5, -2.6]
predict(b_long_cen, pred_dat)$lambda + 
  predict(b_long_cen, pred_dat)$gamma +
  predict(b_long_cen, pred_dat)$alpha * predict(b_long_cen, pred_dat)$mu




# Comparison to bamlss ----------------------------------------------------


bamlss_long <-  bamlss(f_rebamlss, family = "jmFEHLERSUCHE", 
                       data = dat_sigmas$data, timevar = "obstime", 
                       idvar = "id", sampler = FALSE, maxit = 200)
(sp_bamlsl <- bamlss_long$parameters$lambda$s$`s(survtime)`[10])

b_long_bamls <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmas$data, 
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       tau = list("lambda" = list("s(survtime)" = sp_bamlsl)))

plot(bamlss_long, ask = FALSE, model = "lambda")
plot(b_long_bamls, ask = FALSE, model = "lambda")
# Again too smooth
bamlss_long$parameters$alpha
b_long_bamls$parameters$alpha

# True hazard: [-5, -2.6]
predict(bamlss_long, pred_dat)$lambda + 
  predict(bamlss_long, pred_dat)$gamma +
  predict(bamlss_long, pred_dat)$alpha * predict(bamlss_long, pred_dat)$mu


# Recompute on the centered data
bamlss_long_cen <- bamlss(f_rebamlss, family = "jmFEHLERSUCHE", 
                          data = dat_cen$data, timevar = "obstime",
                          idvar = "id", sampler = FALSE, maxit = 200)
bamlss_long_cen$parameters$gamma$p
bamlss_long_cen$parameters$alpha$p
bamlss_long_cen$parameters$mu$p

# True hazard: [-5, -2.6]
# Centering does not change the model prediction
predict(bamlss_long_cen, pred_dat)$lambda + 
  predict(bamlss_long_cen, pred_dat)$gamma +
  predict(bamlss_long_cen, pred_dat)$alpha * predict(bamlss_long_cen, 
                                                     pred_dat)$mu



# Linear BH, Large Sigma --------------------------------------------------


set.seed(1808)

# Linear hazard, no association, constant marker values
dat_sigmal <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                         lambda = function(t, x) -5 + 0.02*t,
                         gamma = function(x) 0,
                         alpha = list(function(t, x) 0*t),
                         mu = list(function(t, x) 1.25),
                         sigma = function(t, x) 0.001 + 0*t,
                         full = TRUE)
ggplot(dat_sigmal$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(data = dat_sigmal$data %>%
               distinct(id, mu, survtime, event),
             aes(x = survtime, y = mu, shape = factor(event)))

# Different model fits
b_longls <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmal$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200,)
b_longls_spl <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmal$data,
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       tau = list("lambda" = list("s(survtime)" = 0.005)))
b_longls_spu <- bamlss(f_re, family = mjm_bamlss, data = dat_sigmal$data,
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       tau = list("lambda" = list("s(survtime)" = 5000)))

# Baseline estimates worsen
# b_longls:  (-2.5, 2) -- 4.5/120 = 0.038
# b_longls_spl: (-1.5, 2.5) -- 4/120 = 0.033
# b_longls_spu: (-4, 4) -- 8/120 = 0.07
plot(b_longls, ask = FALSE, model = "lambda")
plot(b_longls_spl, ask = FALSE, model = "lambda")
plot(b_longls_spu, ask = FALSE, model = "lambda")

# Parameter estimates also worsen
c(b_long$parameters$gamma[[1]] + 
    b_long$parameters$alpha[[1]]*b_long$parameters$mu$p[[1]],
  b_longls$parameters$gamma[[1]] + 
    b_longls$parameters$alpha[[1]]*b_longls$parameters$mu$p[[1]])

# True hazard: [-5, -2.6]
predict(b_longls, pred_dat)$lambda + 
  predict(b_longls, pred_dat)$gamma +
  predict(b_longls, pred_dat)$alpha * predict(b_longls, pred_dat)$mu
predict(b_longls_spl, pred_dat)$lambda + 
  predict(b_longls_spl, pred_dat)$gamma +
  predict(b_longls_spl, pred_dat)$alpha * predict(b_longls_spl, pred_dat)$mu
predict(b_longls_spu, pred_dat)$lambda + 
  predict(b_longls_spu, pred_dat)$gamma +
  predict(b_longls_spu, pred_dat)$alpha * predict(b_longls_spu, pred_dat)$mu



# Constant BH, Small Sigma ------------------------------------------------


set.seed(1808)

# Linear hazard, no association, constant marker values
dat_cbh <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                      lambda = function(t, x) -5,
                      gamma = function(x) 0,
                      alpha = list(function(t, x) 0*t),
                      mu = list(function(t, x) 1.25),
                      sigma = function(t, x) -50 + 0*t,
                      full = TRUE)
ggplot(dat_cbh$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))

# Different model fits
cbh_nolong <- bamlss(f_re, family = mjm_bamlss, data = dat_cbh$data, 
                     timevar = "obstime", sampler = FALSE, maxit = 200,
                     opt_long = FALSE)
cbh_long <- bamlss(f_re, family = mjm_bamlss, data = dat_cbh$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200)
cbh_long$parameters$lambda$s$`s(survtime)`[10]

cbh_long_spl <- bamlss(f_re, family = mjm_bamlss, data = dat_cbh$data,
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       tau = list("lambda" = list("s(survtime)" = 0.005)))
cbh_long_spu <- bamlss(f_re, family = mjm_bamlss, data = dat_cbh$data,
                       timevar = "obstime", sampler = FALSE, maxit = 200,
                       tau = list("lambda" = list("s(survtime)" = 5000)))

plot(cbh_nolong, ask = FALSE, model = "lambda")
plot(cbh_long, ask = FALSE, model = "lambda")
plot(cbh_long_spl, ask = FALSE, model = "lambda")
plot(cbh_long_spu, ask = FALSE, model = "lambda")

c(cbh_nolong$parameters$gamma[[1]], cbh_long$parameters$gamma[[1]])
predict(cbh_nolong, pred_dat)$lambda + 
  predict(cbh_nolong, pred_dat)$gamma
predict(cbh_long, pred_dat)$lambda + 
  predict(cbh_long, pred_dat)$gamma +
  predict(cbh_long, pred_dat)$alpha * predict(cbh_long, pred_dat)$mu



# Linear BH, Small Sigma, Association -------------------------------------


set.seed(1808)

# Linear hazard, constant association, constant marker values
dat_assocs <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                         lambda = function(t, x) - 6.25 + 0.02*t,
                         gamma = function(x) 0,
                         alpha = list(function(t, x) 1 + 0*t),
                         mu = list(function(t, x) 1.25),
                         sigma = function(t, x) -50 + 0*t,
                         full = TRUE)
ggplot(dat_assocs$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_point(aes(x = survtime, shape = factor(event)))


# Different model specifications
a_nolong <- bamlss(f_re, family = mjm_bamlss, data = dat_assocs$data, 
                   timevar = "obstime", sampler = FALSE, maxit = 200, 
                   opt_long = FALSE)
a_long <- bamlss(f_re, family = mjm_bamlss, data = dat_assocs$data, 
                 timevar = "obstime", sampler = FALSE, maxit = 200)
a_long_sml <- bamlss(f_re, family = mjm_bamlss, data = dat_assocs$data, 
                     timevar = "obstime", sampler = FALSE, maxit = 200,
                     tau = list("lambda" = list("s(survtime)" = 0.005)))
a_long_smu <- bamlss(f_re, family = mjm_bamlss, data = dat_assocs$data, 
                     timevar = "obstime", sampler = FALSE, maxit = 200,
                     tau = list("lambda" = list("s(survtime)" = 5000)))
a_bamlss <- bamlss(f_rebamlss, family = "jm", data = dat_assocs$data, 
                   timevar = "obstime", idvar = "id", sampler = FALSE, 
                   maxit = 200)

# BH plots, best approximation to truth with sml
plot(a_nolong, model = "lambda", ask = FALSE)
plot(a_long, model = "lambda", ask = FALSE)
plot(a_long_sml, model = "lambda", ask = FALSE)
plot(a_long_smu, model = "lambda", ask = FALSE)
plot(a_bamlss, model = "lambda", ask = FALSE)

# True hazard: [-5, -2.6]
predict(a_nolong, pred_dat)$lambda + 
  predict(a_nolong, pred_dat)$gamma
predict(a_long, pred_dat)$lambda + 
  predict(a_long, pred_dat)$gamma +
  predict(a_long, pred_dat)$alpha * predict(a_long, pred_dat)$mu
predict(a_long_sml, pred_dat)$lambda + 
  predict(a_long_sml, pred_dat)$gamma +
  predict(a_long_sml, pred_dat)$alpha * predict(a_long_sml, pred_dat)$mu
predict(a_long_smu, pred_dat)$lambda + 
  predict(a_long_smu, pred_dat)$gamma +
  predict(a_long_smu, pred_dat)$alpha * predict(a_long_smu, pred_dat)$mu
predict(a_bamlss, pred_dat)$lambda + 
  predict(a_bamlss, pred_dat)$gamma +
  predict(a_bamlss, pred_dat)$alpha * predict(a_bamlss, pred_dat)$mu

model_list <- c("a_nolong", "a_long", "a_long_sml", "a_long_smu", "a_bamlss")
sapply(model_list, function(x) get(x)$parameters$gamma[[1]])
sapply(model_list, function(x) get(x)$parameters$alpha[[1]])




# Constant Association with Centered Data ---------------------------------


dat_assocc <- dat_assocs
dat_assocc$data$y <- dat_assocc$data$y - mean(dat_assocc$data$y)

a_cen <- bamlss(f_re, family = mjm_bamlss, data = dat_assocc$data, 
                timevar = "obstime", sampler = FALSE, maxit = 200)
a_cen_bamlss <- bamlss(f_rebamlss, family = "jm", data = dat_assocc$data, 
                       timevar = "obstime", idvar = "id", sampler = FALSE, 
                       maxit = 200)

# BH is deeply affected
plot(a_cen, model = "lambda", ask = FALSE)
plot(a_bamlss, model = "lambda", ask = FALSE)

sapply(c("a_cen", "a_cen_bamlss"), function(x) get(x)$parameters$gamma[[1]])
sapply(c("a_cen", "a_cen_bamlss"), function(x) get(x)$parameters$alpha[[1]])



# Constant Association with Linear Mu -------------------------------------


set.seed(1808)

# Linear hazard, constant association, constant marker values
dat_lin <- simMultiJM(nsub = 30, nmark = 1, M = 1,
                      lambda = function(t, x) - 6.25 + 0.015*t,
                      gamma = function(x) 0,
                      alpha = list(function(t, x) 1 + 0*t),
                      mu = list(function(t, x) 1.25 + 0.005*t),
                      sigma = function(t, x) -50 + 0*t,
                      full = TRUE)
ggplot(dat_lin$data, aes(x = obstime, y = y, group = id)) + 
  geom_line() +
  geom_segment(dat_lin$data %>% group_by(id) %>% 
                 filter(obstime == max(obstime)),
               mapping = aes(x = obstime, y = y, xend = survtime, yend = y),
               linetype = "dotted") +
  geom_point(dat_lin$data %>% group_by(id) %>% 
               filter(obstime == max(obstime)),
             mapping = aes(x = survtime, shape = factor(event)))


# Different model specs
a_lin <- bamlss(f_re, family = mjm_bamlss, data = dat_lin$data,
                timevar = "obstime", sampler = FALSE, maxit = 200)
a_lin_sml <- bamlss(f_re, family = mjm_bamlss, data = dat_lin$data, 
                    timevar = "obstime", sampler = FALSE, maxit = 200,
                    tau = list("lambda" = list("s(survtime)" = 0.005)))
a_lin_smu <- bamlss(f_re, family = mjm_bamlss, data = dat_lin$data, 
                    timevar = "obstime", sampler = FALSE, maxit = 200,
                    tau = list("lambda" = list("s(survtime)" = 5000)))
a_lin_bamlss <- bamlss(f_rebamlss, family = "jm", data = dat_lin$data, 
                   timevar = "obstime", idvar = "id", sampler = FALSE, 
                   maxit = 200)

plot(a_lin, model = "lambda", ask = FALSE)
plot(a_lin_sml, model = "lambda", ask = FALSE)
plot(a_lin_smu, model = "lambda", ask = FALSE)
plot(a_lin_bamlss, model = "lambda", ask = FALSE)

# Bamlss estimates alpha parameter to 0
model_list <- c("a_lin", "a_lin_sml", "a_lin_smu", "a_lin_bamlss")
sapply(model_list, function(x) get(x)$parameters$gamma[[1]])
sapply(model_list, function(x) get(x)$parameters$alpha[[1]])
sapply(model_list, function(x) get(x)$parameters$mu[[2]][[1]])
sapply(model_list, function(x) get(x)$parameters$mu[[2]][[2]])


# True hazard: (-5, -2.6)
predict(a_lin, pred_dat)$lambda + 
  predict(a_lin, pred_dat)$gamma +
  predict(a_lin, pred_dat)$alpha * predict(a_lin, pred_dat)$mu
predict(a_lin_sml, pred_dat)$lambda + 
  predict(a_lin_sml, pred_dat)$gamma +
  predict(a_lin_sml, pred_dat)$alpha * predict(a_lin_sml, pred_dat)$mu
predict(a_lin_smu, pred_dat)$lambda + 
  predict(a_lin_smu, pred_dat)$gamma +
  predict(a_lin_smu, pred_dat)$alpha * predict(a_lin_smu, pred_dat)$mu
predict(a_lin_bamlss, pred_dat)$lambda + 
  predict(a_lin_bamlss, pred_dat)$gamma +
  predict(a_lin_bamlss, pred_dat)$alpha * predict(a_lin_bamlss, pred_dat)$mu
# bamlss has the smallest spread



# Use PCRE Model Terms ----------------------------------------------------

mfpcaT <- create_true_MFPCA(M = 1, nmarker = 1)
f_pcre <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ 1 + obstime + s(id, wfpc, bs = "pcre", xt = list("mfpc" = mfpcaT)),
  sigma ~ 1,
  alpha ~ 1
)

# Different model specs
pcre_a_lin <- bamlss(f_pcre, family = mjm_bamlss, data = dat_lin$data,
                     timevar = "obstime", sampler = FALSE, maxit = 200)
pcre_a_lin_sml <- bamlss(f_pcre, family = mjm_bamlss, data = dat_lin$data,
                         timevar = "obstime", sampler = FALSE, maxit = 200,
                         tau = list("lambda" = list("s(survtime)" = 0.005)))
pcre_a_lin_smu <- bamlss(f_pcre, family = mjm_bamlss, data = dat_lin$data, 
                         timevar = "obstime", sampler = FALSE, maxit = 200,
                         tau = list("lambda" = list("s(survtime)" = 5000)))

# Only small differences to the a_lin models based on RE
plot(pcre_a_lin, model = "lambda", ask = FALSE)
plot(pcre_a_lin_sml, model = "lambda", ask = FALSE)
plot(pcre_a_lin_smu, model = "lambda", ask = FALSE)

