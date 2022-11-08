# Also folgendes:
# JM_example führt zu numerischen Problemen. Interessanterweise scheint das aber
# an der Datensimulation zu liegen. Auch den "Simple" Fall unten, wo nur eine 
# Dimension generiert wird, kann man nicht mit bamlss berechnen. Da muss also 
# noch irgendwo ein Problem in SimMultiJM liegen. Deshalb versuche jetzt ähnlich
# zu Meikes Simulation einen Datensatz mit MJM zu erstellen (univariates JM),
# damit man sicherstellen kann, dass die Daten passen sollten.


# Meike's Simulation

# data_setting "a"
# long_setting "linear"
# alpha_setting "constant"
# dalpha_setting "zero"
# nsub 150
# nreps 1

nsub <- 150
times <- seq(0, 120, 1)
probmiss <- 0.75
long_setting <- "linear"
alpha_setting <- "constant"
dalpha_setting <- "zero"
seed <- 2505 + 1


dat_m <- bamlss::simJM(nsub = nsub, times = times, probmiss = probmiss,
                       long_setting = long_setting,
                       alpha_setting = alpha_setting, 
                       dalpha_setting = dalpha_setting,
                       seed = seed, full = TRUE)

# The following covariates are generated:
# x1 <- runif(nsub, -3, 3) 
# x2 <- runif(nsub, -3, 3) 
# x3 <- rbinom(nsub, 1, 0.5)

# The following random effects are generated:
# r1 <- rnorm(nsub, 0, 0.25) 
# r2 <- rnorm(nsub, 0, 0.4) 




# Own implementation
source("R/simMultiJM.R")

# Das hier wäre das Ziel
dat_a <- simMultiJM(nsub = nsub, times = times, probmiss = probmiss, 
                    maxfac = 1.5, nmark = 1, M = 2, ncovar = 2, 
                    lambda = function (time, x) {
                      1.4 * log((time + 10) / 1000) - 1.5
                    },
                    gamma = function(x) {
                      0.3*x[, 1]
                    },
                    alpha = list(function(time, x) {
                      1 + 0*time
                    }),
                    mu = list(function (time, x) {
                      1.25 + 0.6*sin(x[, 2]) + (-0.01)*time# + 
                        #r[, 1] + r[, 2]*0.02*time
                    }),
                    sigma = function(time, x) {
                      log(0.3) + 0*time
                    }, tmax = NULL, seed = 1808,
                    mfpc_args = list(type = "split", eFunType = "Poly",
                                     ignoreDeg = NULL, eValType = "linear",
                                     eValScale = 1),
                    full = FALSE, file = NULL)

# Was machn wa mit Random Intercept / Random Slope?
# Erst ma ignoriern ma des

# Kann man das Modell berechnen?
f_mjm <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1 + x1,
  mu ~ obstime + s(x2) + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
f_jm <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1 + x1,
  mu ~ obstime + s(x2) + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
library(bamlss)
# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")
source("R/survint.R")
b_jm <- bamlss(f_jm , family = "jm", data = dat_a, timevar = "obstime",
               idvar = "id")
sink("b_simdat.txt")
b_mjm <- bamlss(f_mjm, family = mjm_bamlss, data = dat_a, timevar = "obstime")
sink()


# Weißt was, jetzt fangen wir erst mal an mit den Basics:
# Was passiert eigentlich in
d_simp <- simMultiJM(nsub = 300, times = seq(0, 1, length.out = 121),
                     probmiss = 0.75, maxfac = 1.5, nmark = 1, M = 6, 
                     ncovar = 2,
                     lambda = function(t, x) {
                       1.4*log((120*t + 10)/1000)
                     },
                     alpha = list(function(t, x) {
                       0.3 + 0*t
                     }),
                     mu = list(function(t, x){
                       1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
                     }),
                     sigma = function(t, x) {-50 + 0*t}, 
                     full = TRUE,
                     mfpc_args = list(type = "split", eFunType = "Poly",
                                      ignoreDeg = NULL, eValType = "linear",
                                      eValScale = 1))
