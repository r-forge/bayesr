
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


# Center Longitudinal Data ------------------------------------------------

dat_cen <- dat_sigmas
dat_cen$data$y <- dat_cen$data$y - mean(dat_cen$data$y)

b_long_cen <- bamlss(f_re, family = mjm_bamlss, data = dat_cen$data, 
                     timevar = "obstime", sampler = FALSE, maxit = 200,)
b_long_cen$parameters$gamma$p
b_long_cen$parameters$alpha$p
b_long_cen$parameters$mu$p
b_long_cen$parameters$gamma[[1]] + 
  b_long_cen$parameters$alpha[[1]]*b_long_cen$parameters$mu$p[[1]]


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



# Multivariate Example ----------------------------------------------------

set.seed(1808)
d <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121), 
                lambda = function(t, x) {
                  1.4*log((120*t + 10)/1000)
                })

f_re_mul <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + ti(obstime, by = marker) +
    ti(id, bs = "re") +
    ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b <- bamlss(f_re_mul, family = mjm_bamlss, data = d, timevar = "obstime",
            sampler = FALSE, maxit = 50)


f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_f <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime",
              sampler = FALSE)

MFPCA <- create_true_MFPCA(M = 2, nmarker = 2)
plot(b_f, model = "lambda", ask = FALSE)
