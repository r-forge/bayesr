library(bamlss)
library(funData)
library(tidyverse)
source("R/simJM.R")
sim_dat <- simJM(nsub = 10, long_setting = "fpc", alpha_setting = "constant",
                 long_df = 3, sigma = 0.05, seed = 120, eval_scale = 2)

# The data are generated so that in Yao (2007) notation
# Formula 6: h_i(t) = h_0(t)*exp(\gamma X_i(t) + V_i(t)'\zeta)
# 
# with
# t = sim_dat$obstime
# h_0(t) = 1.4*log((time + 10)/1000) - 1.5 ---- function lambda()
# \gamma = 1 ---- function alpha()
# V_i(t) = sim_dat$x1 ---- function gamma()
# \zeta = 0.3 ---- function gamma()
# 
# Formula 1: X_i(t_ij) = \mu(t_ij) + Z_i(t_ij)'\beta + 
#               \sum_{k=1}^K \xi_ik\phi_k(t_ij) + \epsilon_ij ---- function mu()
# \mu_ij = 1.25
# \Z_i(t_ij)' = (sin(sim_dat$x2), obstime)
# \beta = (0.6, -0.01)
# --------------------- NOTE -------------------
# The covariate effect of x2 is actually 0.6*sin(x2). This is important as the
# proposed method of Yao is not able to estimate the nonlinear part, only the
# coefficient. It assumes that the transformation function of x2 is fully known.
# ----------------------------------------------
# K = 3
# \xi_ik ~ N(0, \lambda_k) with 
# \lambda_1 = 2*1, \lambda_2 = 2*2/3, \lambda_3 = 2*1/3 ---- function gen_fpc()
# ---- eVal()
# \phi_1, \phi_2, \phi_3 --- function mu() ---- eFun() ---- funData:::efPoly()
# \epsilon_ij ~ N(0, 0.05)
 
# Short data frame for events
sim_dat_short <- sim_dat %>%
  group_by(id) %>%
  slice(tail(row_number(), 1)) %>%
  ungroup()

# Plot the data
ggplot(data = sim_dat, aes(x = obstime, y = y, color = id)) + 
  geom_line() +
  geom_point(data = sim_dat_short, 
             mapping = aes(x = survtime, y = y, shape = factor(event)),
             size = 2) +
  scale_shape_manual(values = c(1, 4))



# Simulate a larger data set ----------------------------------------------


source("R/simMultiJM.R")
# Generate data with independent random intercepts
dat_uni <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                      probmiss = 0.75, maxfac = 2,
                      nmark = 1, param_assoc = FALSE, M = 4, 
                      FPC_bases = NULL, FPC_evals = NULL, 
                      mfpc_args = list("type" = "split",
                                       "eFunType" = "PolyHigh",
                                       "ignoreDeg" = 1, 
                                       "eValType" = "linear",
                                       "eValScale" = 3),
                      re_cov_mat = NULL, ncovar = 2,
                      lambda = function(t, x) {
                        1.27 * t^(0.27)
                      },
                      gamma = function(x) {
                        - 7.8 + 0.48*x[, 3]
                      },
                      alpha = list(function(t, x) {
                        0.64 + 0*t
                      }),
                      mu = list(function(t, x, r){
                        2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                      }),
                      sigma = function(t, x) {
                        log(0.6) + 0*t
                      }, 
                      tmax = NULL, seed = 1808, 
                      full = TRUE, file = NULL)
table(dat_uni$data_short$event)
ggplot(data = dat_uni$data_short %>% mutate(row = as.numeric(row.names(.)),
                                            event = factor(event)), 
       aes(x = row, y = survtime, shape = event, color = event)) +
  geom_point()

ggplot(dat_uni$data, aes(x = obstime, y = y, color = id)) +
  geom_point() +
  geom_line(aes(y = mu)) +
  theme(legend.position = "none")
sim_data <- dat_uni$data
save(sim_data, file = "inst/objects/sim_uni_fpc.Rdata")
