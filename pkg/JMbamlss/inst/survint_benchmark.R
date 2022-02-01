

# Bivariate Joint Model Example For Optimizing Calculations ---------------


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
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")


# Data Generation ---------------------------------------------------------

# if(!exists("d")) {
#   
#   library("survival")
#   library("bamlss")
#   library("MFPCA")
#   
#   
#   set.seed(1808)
#   d <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121),
#                   lambda = function(t, x) {
#                     1.4*log((120*t + 10)/1000)
#                   },
#                   alpha = rep(list(function(t, x) {
#                     0.3 + 0*t
#                   }), 2),
#                   mu = rep(list(function(t, x){
#                     1.25 + 0.6*sin(x[, 2]) + (-0.01)*t
#                   }), 2),
#                   sigma = function(t, x) {-50 + 0*t})
#   
# }
# 
# f <- list(
#   Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
#   gamma ~ 1,
#   mu ~ s(obstime) + s(id, bs = "re"),
#   sigma ~ 1,
#   alpha ~ 1
# )
# 
# b <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime")
# save(update_obj, file = "inst/objects/update_obj.Rdata")

# Take values from last iteration (code above)
load("inst/objects/update_obj.Rdata")
for ( i in names(update_obj)){
  assign(i, update_obj[[i]])
}

# Try Out Updating Objects ------------------------------------------------


# Lambda ------------------------------------------------------------------

# Updating lambda
newstate <- update_mjm_lambda(x$lambda$smooth.construct[[1]], y = y, nu = nu,
                              eta = eta, eta_timegrid = eta_timegrid,
                              survtime = survtime)

# Survival interval
int_i <- survint_gq(pred = "lambda", pre_fac = exp(eta$gamma),
                    omega = exp(eta_timegrid),
                    int_vec = x$lambda$smooth.construct[[1]]$Xgrid, 
                    weights = attr(y, "gq_weights"),
                    survtime = survtime)


# Gamma -------------------------------------------------------------------

# Updating gamma
newstate <- update_mjm_gamma(x$gamma$smooth.construct[[1]], y = y, nu = nu,
                             eta = eta, eta_timegrid = eta_timegrid,
                             survtime = survtime)

# Survival interval
int_i <- survint_gq(pred = "gamma", pre_fac = exp(eta$gamma), 
                    pre_vec = x$gamma$smooth.construct[[1]]$X,
                    omega = exp(eta_timegrid),
                    weights = attr(y, "gq_weights"),
                    survtime = survtime)


# Alpha -------------------------------------------------------------------

# Updating alpha
newstate <- update_mjm_alpha(x$alpha$smooth.construct[[1]], y = y, 
                             nu = nu,
                             eta = eta, eta_timegrid = eta_timegrid,
                             eta_timegrid_mu = eta_timegrid_mu,
                             eta_T_mu = eta_T_mu, survtime = survtime)

# Survival interval
int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                    omega = exp(eta_timegrid),
                    int_fac = eta_timegrid_mu, 
                    int_vec = x$alpha$smooth.construct[[1]]$Xgrid,
                    weights = attr(y, "gq_weights"),
                    survtime = survtime)



# Mu ----------------------------------------------------------------------

# Updating mu
newstate <- update_mjm_mu(x$mu$smooth.construct[[1]], y = y, nu = nu,
                          eta = eta, eta_timegrid = eta_timegrid,
                          eta_timegrid_alpha = eta_timegrid_alpha,
                          survtime = survtime)

# Survival interval
int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                    omega = exp(eta_timegrid),
                    int_fac = eta_timegrid_alpha, 
                    int_vec = x$mu$smooth.construct[[1]]$Xgrid,
                    weights = attr(y, "gq_weights"),
                    survtime = survtime)

# Random effects
int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                    omega = exp(eta_timegrid),
                    int_fac = eta_timegrid_alpha, 
                    int_vec = x$mu$smooth.construct[[2]]$Xgrid,
                    weights = attr(y, "gq_weights"),
                    survtime = survtime)
