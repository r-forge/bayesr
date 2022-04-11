
# How to plot fitted REs --------------------------------------------------

# NOTE: new version of bamlss needed!
library(survival)
library(bamlss)
library(MFPCA)
library(tidyverse)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

# Load simulated data and fitted model
load("inst/objects/find_sim_data.Rdata")
load("inst/objects/find_sim.Rdata")

# Construct true underlying FPCs
mfpca <- create_true_MFPCA(M = 8, nmarker = 2, argvals = seq(0, 25, by = 0.25),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(80:73))


# Data for prediction -----------------------------------------------------

# Which ids to plot
set.seed(188)
ids <- sample(1:150, 5)

# Set up data for prediction, start with one dimension
times <- seq(0, 20, 1)
newdat <- d_indepri$data_short[rep(1:150, each = length(times)), ]
newdat$obstime <- rep(times, 150)
newdat <- newdat[newdat$survtime > newdat$obstime, ]

# Expand to two dimensions
wfpcs <- eval_mfpc(mfpca, timepoints = newdat$obstime)
newdat <- rbind(newdat,
                newdat %>% mutate(marker = fct_recode(marker, "m2" = "m1")))
newdat <- cbind(newdat, wfpcs) %>% filter(id %in% ids)



# Prediction --------------------------------------------------------------

# Predict functional random intercepts with predict
newdat$fri_pred <- predict(b_sim, newdat, model = "mu", term = "id", 
                           intercept = FALSE)

# Predict functional random intercepts manually
newdat_id <- newdat %>% droplevels() %>% split(.$id)
pars_id <- lapply(ids, function(i) (1+(i-1)*8):(1+(i-1)*8+7))
samp_means <- summary(b_sim$samples, quantiles = 0.5)$quantiles
newdat_id <- do.call(rbind, mapply(function(dat, pars) {
  ps <- samp_means[12:1211][pars]
  dat$fri_man <- c(as.matrix(dat[, 25:32]) %*% ps)
  dat
}, dat = newdat_id, pars = pars_id, SIMPLIFY = FALSE))



# Plot predictions --------------------------------------------------------

# Plot functional random effects
ggplot(data = newdat_id, aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fri_pred), size = 1.2) +
  geom_line(aes(y = fri_man), size = 1.2, linetype = "dashed") +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)
# With predict function they are the same over the two markers

