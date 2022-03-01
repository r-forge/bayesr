library(Matrix)
library(mvtnorm)
library("survival")
library("bamlss")
library("MFPCA")

# Alternative bamlss Code (more options)
source("Fehlersuche/JM.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")
source("R/survint.R")

load("inst/objects/simu_meike/a_linear_constant_zero_150_1.RData")


# Optimize in BAMLSS -----------------------------------------------------

seed <- 5555 + 1
set.seed(seed)

f1 <- list(                                                                  
  Surv2(survtime, event, obs = y) ~ s(survtime, k = 6, bs="ps"),
  gamma ~ s(x1, k = 6, bs="ps"),
  mu ~ obstime + s(x2, k = 8, bs="ps") + s(id, bs="re"),
  sigma ~ 1,
  alpha ~ 1, 
  dalpha ~ -1
)

# Fit posterior mode
a <- tryCatch({bamlss(f1, data = d$data, family = "jm", timevar = "obstime",
                      idvar = "id", subdivisions = 25, 
                      maxit = 1500, sampler=FALSE, do.optim2 = FALSE)},
              error = function(cond){
                message(cond)
                return(2)})


# Sampler in BAMLSS -------------------------------------------------------

# b1 <- bamlss(f1, data = d$data, family = "jmFEHLERSUCHE", timevar = "obstime",
#              idvar = "id", subdivisions = 25,# n.iter = 20, burnin = 2,
#              prop_pred = "lambda", optimizer = FALSE, start = parameters(a),
#              verbose_sampler = TRUE, verbose = FALSE)
sink("bamlss.txt")
b0 <- bamlss(f1, data = d$data, family = "jmFEHLERSUCHE", timevar = "obstime",
             idvar = "id", subdivisions = 250, n.iter = 20, burnin = 0,
             #prop_pred = "lambda",
             optimizer = FALSE, start = parameters(a),
             verbose_sampler = TRUE, verbose = FALSE, prop_list = TRUE)
sink()


# Fit model in MJM --------------------------------------------------------

seed <- 5555 + 1
set.seed(seed)

f2 <- list(                                                                  
  Surv2(survtime, event, obs = y) ~ s(survtime, k = 6, bs="ps"),
  gamma ~ s(x1, k = 6, bs="ps"),
  mu ~ obstime + s(x2, k = 8, bs="ps") + s(id, bs="re"),
  sigma ~ 1,
  alpha ~ 1
)

# Fit posterior mode
sink("MJM_prop.txt")
b2 <- bamlss(f2, data = d$data, family = mjm_bamlss, timevar = "obstime", 
             n.iter = 20, burnin = 0, subdivisions = 200,
             #prop_pred = "lambda", 
             optimizer = FALSE, start = parameters(a), verbose_sampler = TRUE,
             prop_list = prop_list)
sink()
