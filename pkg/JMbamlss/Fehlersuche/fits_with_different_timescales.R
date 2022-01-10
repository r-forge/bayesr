
# Comparison between data examples ----------------------------------------

# Data Generation
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

# Load Packages
library("refund")
library("bamlss")
library("MFPCA")

Surv2 <- bamlss:::Surv2

par(mfrow = c(1,2))

set.seed(1808)


d_new_full <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121), 
                         lambda = function(t, x) {
                           1.4*log((120*t + 10)/1000)
                          }, M = 2, full = TRUE)
d_new <- d_new_full$data
d_new120 <- d_new
d_new120$survtime <- d_new120$survtime*120
d_new120$obstime <- d_new120$obstime*120

d_old_full <- simMultiJM(nsub = 50,
                         lambda = function(t, x) {
                           1.4*log((t + 10)/1000)
                         }, M = 2, full = TRUE)
d_old <- d_old_full$data
d_old1 <- d_old
d_old1$survtime <- d_old1$survtime/120
d_old1$obstime <- d_old1$obstime/120

# True functional principal components
mfpca1 <- create_true_MFPCA(M = 2, nmarker = 2, 
                            argvals = seq(0, 1, length.out = 121))
mfpca120 <- create_true_MFPCA(M = 2, nmarker = 2)

# The baseline hazard is centered so the mean of the observed eta_lambdas
# is added to the gamma predictor
mean_new <- mean(attr(d_new_full$data_short, "f_lambda"))
mean_old <- mean(attr(d_old_full$data_short, "f_lambda"))


# PCRE Model for new example ----------------------------------------------


f_new <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ x1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = mfpca1)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_new <- bamlss(f_new, family = mjm_bamlss, data = d_new, timevar = "obstime",
                sampler = FALSE)


curve(1.4*log((120*x + 10)/1000) - mean_new, from = 0, to = 1)
plot(b_new, model = "lambda", ask = FALSE)



# Change the time-scale ---------------------------------------------------


f_new120 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ x1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = mfpca120)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_new120 <- bamlss(f_new120, family = mjm_bamlss, data = d_new120, 
                   timevar = "obstime", sampler = FALSE)


curve(1.4*log((x + 10)/1000) - mean_new, from = 0, to = 120)
plot(b_new120, model = "lambda", ask = FALSE)



# Old example -------------------------------------------------------------


f_old <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ x1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = mfpca120)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

# Model fit with PCRE and plot
b_old <- bamlss(f_old, family = mjm_bamlss, data = d_old, timevar = "obstime",
                sampler = FALSE)

curve(1.4*log((x + 10)/1000) - mean_old, from = 0, to = 120)
plot(b_old, model = "lambda", ask = FALSE)



# Old data set with FREs on 120 -------------------------------------------


# Model fit with FRE and plot
f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ x1,
  mu ~ -1 + marker + ti(obstime, by = marker) +
    ti(id, bs = "re") +
    ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_re <- bamlss(f_re, family = mjm_bamlss, data = d_old, timevar = "obstime",
               sampler = FALSE)
curve(1.4*log((x + 10)/1000) - mean_old, from = 0, to = 120)
plot(b_re, model = "lambda", ask = FALSE)




# Old data set with first reducing time scale -----------------------------

# List element has to be updated
f_old1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = mfpca1)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

# Model fit with pcre and plot
b_old1 <- bamlss(f_old1, family = mjm_bamlss, data = d_old1, 
                 timevar = "obstime", sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000) - mean_old, from = 0, to = 1)
plot(b_old1, model = "lambda", ask = FALSE)



# Old data set with FREs on 1 ---------------------------------------------


# Model fit with fre and plot
b_re1 <- bamlss(f_re, family = mjm_bamlss, data = d_old1, timevar = "obstime",
                 sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000) - mean_old, from = 0, to = 1)
plot(b_re, model = "lambda", ask = FALSE)
