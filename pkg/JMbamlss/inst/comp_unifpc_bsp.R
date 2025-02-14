
# Set up ------------------------------------------------------------------

library(bamlss)
library(refund)
library(funData)
library(tidyverse)
source("R/simJM.R")
sim_dat <- simJM(nsub = 100, long_setting = "fpc", alpha_setting = "constant",
                 long_df = 3, sigma = 0.05, seed = 120, eval_scale = 2)

# Simulated data follow the model
# lambda = 1.4*log((obstime + 10)/1000) - 1.5
# gamma = 0.3*x1
# mu = 1.25 + 0.6*sin(x2) + (-0.01)*obstime + 3fpcs_poly_linear_scale2
# sigma = 0.05
# alpha = 1

# Plot the data

sim_dat_short <- sim_dat %>%
  group_by(id) %>%
  slice(tail(row_number(), 1)) %>%
  ungroup()
ggplot(data = sim_dat, aes(x = obstime, y = y, color = id)) + 
  geom_line() +
  geom_point(data = sim_dat_short, 
             mapping = aes(x = survtime, y = y, shape = factor(event)),
             size = 2) +
  scale_shape_manual(values = c(1, 4)) +
  guides(color = FALSE)


# JM based on univariate FPCs ---------------------------------------------

f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ x1,
  mu ~ obstime + s(x2, bs = "ps") + s(id, wfpc.1, wfpc.2, wfpc.3, bs = "pcre"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
# the following model can also be loaded from inst/objects
m <- bamlss(f, data = sim_dat, family = "jm", timevar = "obstime", idvar = "id")
# save(m, file = "objects/m_fpc.Rdata")
# optimizer ca. 2 min
# sampler ca. 20 min
nrow(coef(m))
# 331 parameters



# Alternative model based on B-spline RE ----------------------------------

f1 <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ x1,
  mu ~ obstime + s(x2, bs = "ps") + ti(id, bs = "re") +
    ti(id, obstime, bs = c("re", "cr"), k = c(nlevels(sim_dat$id), 8)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
# the following model can also be loaded from inst/objects
m1 <- bamlss(f1, data = sim_dat, family = "jm", timevar = "obstime",
             idvar = "id")
# save(m1, file = "objects/m_bsp.Rdata")
# optimizer ca 1 min
# sampler ca 8 min
nrow(coef(m1))
# 837 parameters


# Time difference is the other way round when using mgcv
system.time(
  m_gcv <- bam(y ~  obstime + s(x2, bs = "ps") + 
                 s(id, wfpc.1, wfpc.2, wfpc.3, bs = "pcre"),
               data = sim_dat)
)
system.time(
  m1_gcv <- bam(y ~ obstime + s(x2, bs = "ps") + ti(id, bs = "re") +
                  ti(id, obstime, bs = c("re", "cr"), 
                     k = c(nlevels(sim_dat$id), 8)),
                data = sim_dat)
)

# Compare model fits ------------------------------------------------------

# Understand bamlss predict
newdat <- sim_dat[rep(1, 5), ]
newdat$obstime <- 0:4
m1_int <- predict(m1, newdat, term = "1")$mu
m1_obs <- predict(m1, newdat, term = "obstime", intercept = FALSE)$mu
# Why not 0 for obstime=0?
# Answer: obstime is contained in two model terms so it automatically returns
#         the sum of all the model terms containing obstime
m1_x2 <- predict(m1, newdat, term = "s(x2)", intercept = FALSE)$mu
m1_id <- predict(m1, newdat, term = "ti(id)", intercept = FALSE)$mu
m1_fre <- predict(m1, newdat, term = "ti(id,obstime)", intercept = FALSE)$mu
m1_mu <- predict(m1, newdat)$mu

# Only adds up if you leave out functional random effects prediction
all.equal(m1_mu, m1_int+m1_obs+m1_x2+m1_id)#+m1_fre

library(mgcv)

