
# Compare Runtime ---------------------------------------------------------



# SetUp for PBC Example ---------------------------------------------------

library(survival)
library(bamlss)
data("pbc2", package = "JMbayes")

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


## Set up the model formula including
## functional random intercepts using ti().
f_jm <- list(
  Surv2(years, status2, obs = log(serBilir)) ~ s(years),
  gamma ~ s(age) + drug + sex,
  mu ~ ti(id,bs="re") + 
    ti(year),
  sigma ~ 1,
  alpha ~ s(years),
  dalpha ~ -1
)

f_mjm <- list(
  Surv2(years, status2, obs = log(serBilir)) ~ s(years),
  gamma ~ s(age) + drug + sex,
  mu ~ ti(year) + ti(id,bs="re"),
  sigma ~ 1,
  alpha ~ s(years),
  dalpha ~ -1
)



# Runtime Transform Function ----------------------------------------------

system.time(b_jm <- bamlss(f_jm, data = pbc2, family = "jm", timevar = "year", 
                           idvar = "id", 
                           optimizer = FALSE, sampler = FALSE))
system.time(b_mjm <- bamlss(f_mjm, family = mjm_bamlss, data = pbc2, 
                            timevar = "year", 
                            optimizer = FALSE, sampler = FALSE))

times <- microbenchmark(
  bamlss(f_jm, data = pbc2, family = "jm", timevar = "year", 
         idvar = "id", 
         optimizer = FALSE, sampler = FALSE),
  bamlss(f_mjm, family = mjm_bamlss, data = pbc2, 
         timevar = "year", 
         optimizer = FALSE, sampler = FALSE)
)
boxplot(times)



`# Runtime Transform + Update ----------------------------------------------

system.time(b_jm <- bamlss(f_jm, data = pbc2, family = "jm", timevar = "year", 
                           idvar = "id", 
                           maxit = 1, sampler = FALSE))
system.time(b_mjm <- bamlss(f_mjm, family = mjm_bamlss, data = pbc2, 
                            timevar = "year", 
                            maxit = 1, sampler = FALSE))



# Full Model --------------------------------------------------------------

# ATTENTION: THIS RUNS WAY TOO LONG
## Set the seed for reproducibility.
set.seed(123)

## Estimate model.
b_jm <- bamlss(f_jm, data = pbc2, family = "jm", timevar = "year", idvar = "id")
b_mjm <- bamlss(f_mjm, family = mjm_bamlss, data = pbc2, timevar = "year")
            