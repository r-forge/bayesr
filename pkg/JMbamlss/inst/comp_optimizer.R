library(bamlss)
filepath <- "data"

# Generate Data -----------------------------------------------------------

seed <- 2505

nsub <- 150
times <- seq(0, 120, 1)
probmiss <- 0.75

d <- simJM(nsub = nsub, times = times, probmiss = probmiss,
           long_setting = "linear",
           alpha_setting = "constant", 
           dalpha_setting = "zero",
           seed = seed+1, full = TRUE)

name <- paste("a", "linear", "constant", "zero", 150, 1, sep="_")
save(d, file = paste0(filepath, name ,".RData"))



# Fit Model in BAMLSS -----------------------------------------------------

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


              

# Fit Model in MJM --------------------------------------------------------

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")

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
b <- tryCatch({bamlss(f2, data = d$data, family = mjm_bamlss,
                      timevar = "obstime", idvar = "id", subdivisions = 25, 
                      maxit = 1500, sampler = FALSE)},
              error = function(cond){
                message(cond)
                return(2)})
