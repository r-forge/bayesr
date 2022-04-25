
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
source("R/MJM_predict.R")

# Load simulated data and fitted model
#load("inst/objects/find_sim_data.Rdata") MUSS MAN NOCH ANPASSEN
load("inst/objects/find_sim.Rdata")

# Construct true underlying FPCs
mfpca <- create_true_MFPCA(M = 8, nmarker = 2, argvals = seq(0, 25, by = 0.25),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(80:73))


# Data for prediction -----------------------------------------------------

# # Which ids to plot
# set.seed(188)
# ids <- sample(1:150, 5)

# # Set up data for prediction, start with one dimension
# times <- seq(0, 20, 1)
# newdat <- d_indepri$data_short[rep(1:150, each = length(times)), ]
# newdat$obstime <- rep(times, 150)
# newdat <- newdat[newdat$survtime > newdat$obstime, ]
# 
# # Expand to two dimensions
# wfpcs <- eval_mfpc(mfpca, timepoints = newdat$obstime)
# newdat <- rbind(newdat,
#                 newdat %>% mutate(marker = fct_recode(marker, "m2" = "m1")))
# newdat <- cbind(newdat, wfpcs) #%>% filter(id %in% ids)
newdat <- d_indepri$data_full


# Prediction --------------------------------------------------------------

# Predict functional random intercepts with predict
#newdat$fri_pred <- predict(b_sim, newdat, model = "mu", term = "id")
#newdat$mu_pred <- predict(b_sim, newdat, model = "mu")
attr(b_sim$y, "marker_name") <- "marker"
newdat$fri_pred_by <- MJM_predict(b_sim, newdata = newdat, model = "mu",
                                  term = "id")

#newdat$mu_pred_by <- MJM_predict(b_sim, newdata = newdat, model = "mu")

# Predict functional random intercepts manually
newdat_id <- newdat %>% #droplevels() %>% 
  split(.$id)
pars_id <- lapply(1:150, function(i) (1+(i-1)*8):(1+(i-1)*8+7))
  #lapply(ids, function(i) (1+(i-1)*8):(1+(i-1)*8+7))
samp_means <- summary(b_sim$samples, quantiles = 0.5)$quantiles
newdat_id <- do.call(rbind, mapply(function(dat, pars) {
  ps <- samp_means[12:1211][pars]
  dat$fri_man <- c(as.matrix(dat[, 25:32]) %*% ps)
  dat
}, dat = newdat_id, pars = pars_id, SIMPLIFY = FALSE))

# True fRI
newdat_id <- newdat_id %>% 
  mutate(fri_tru = s1*wfpc.1/sqrt(mfpca$values[1]) +
           s2*wfpc.2/sqrt(mfpca$values[2]) + 
           s3*wfpc.3/sqrt(mfpca$values[3]) +
           s4*wfpc.4/sqrt(mfpca$values[4]) +
           s5*wfpc.5/sqrt(mfpca$values[5]) +
           s6*wfpc.6/sqrt(mfpca$values[6]) +
           s7*wfpc.7/sqrt(mfpca$values[7]) + 
           s8*wfpc.8/sqrt(mfpca$values[8]))

newdat_id_new <- d_indepri$data_full %>% 
  mutate(fri_tru = s1*fpc.1 +
           s2*fpc.2 + 
           s3*fpc.3 +
           s4*fpc.4 +
           s5*fpc.5 +
           s6*fpc.6 +
           s7*fpc.7 + 
           s8*fpc.8)
newdat_id_new[1:6053, 41] <- (-1)*newdat_id_new[1:6053, 41]

# Plot predictions --------------------------------------------------------

# Plot functional random effects
ggplot(data = newdat_id, aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fri_tru), linetype = "dotted") +
  geom_line(aes(y = fri_pred_by), size = 1.2, linetype = "dashed") +
  geom_line(aes(y = fri_man)) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)
# With predict function they are the same over the two markers

# Aus irgendeinem Grund bekomm ich immer andere Plots geliefert. Wo ist das
# Problem? Und warum überdecken sich die Predictions manchmal, und manchmal
# auch nicht? Das ist alles noch etwas komisch...

ggplot(d_indepri$data %>% filter(id %in% ids), aes(x = obstime, y = mu,
                                                   colour = id)) +
  geom_line(linetype = "dashed") +
  geom_point(aes(y = y)) +
  facet_grid(~marker) +
  theme(legend.position = "none")


ggplot(data = newdat_id_new, aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fri_tru), linetype = "dotted") +
  # geom_line(data = newdat_id_new, aes(x = obstime, y = fri_tru, colour = id),
  #           linetype = "dashed") +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)
