library(survival)
library(bamlss)
library(MFPCA)
library(tidyverse)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")


# Generate data with independent random intercepts
d_indepri <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                        probmiss = 0.75, maxfac = 2,
                        nmark = 2, param_assoc = FALSE, M = 8, 
                        FPC_bases = NULL, 
                        FPC_evals = c(80:73),
                        mfpc_args = list(type = "split", eFunType = "PolyHigh",
                                         ignoreDeg = 1, eValType = "linear",
                                         eValScale = 1),
                        ncovar = 2,
                        lambda = function(t, x) {
                          1.6 * t^(0.6)
                        },
                        gamma = function(x) {
                          - 10 + 0.48*x[, 3]
                        },
                        alpha = list(function(t, x) {
                          0.64 + 0*t
                        }, function(t, x) {
                          -0.64 + 0*t
                        }),
                        mu = list(function(t, x, r){
                          2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                        }, function(t, x, r){
                          2 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3]
                        }),
                        sigma = function(t, x) {
                          log(0.6) + 0*t
                        }, 
                        tmax = NULL, seed = 1808, 
                        full = TRUE, file = NULL)
save(d_indepri, file = "inst/objects/find_sim_data.Rdata")
plot(d_indepri$data_short$survtime)
ggplot(d_indepri$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")



f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + 
    s(id, wfpc.1, wfpc.2, wfpc.3, wfpc.4, wfpc.5, wfpc.6, wfpc.7, wfpc.8,
      bs = "unc_pcre", xt = list("mfpc" = mfpca)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

mfpca <- create_true_MFPCA(M = 8, nmarker = 2, argvals = seq(0, 25, by = 0.25),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(80:73))

# Helperfunction PCRE
source("R/pcre_smooth.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/compile.R")
compile_alex()

sink("find_sim200422.txt")
b_sim1 <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                timevar = "obstime", maxit = 1200, verbose_sampler = TRUE)
sink()
save(b_sim1, file = "inst/objects/find_sim200422.Rdata")


# Why is sigma not accepted? ----------------------------------------------


b <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
            timevar = "obstime", optimizer = FALSE, start = parameters(b_sim1))
mfpca <- create_true_MFPCA(M = 3, nmarker = 2, argvals = seq(0, 25, by = 0.25),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(80:78))
f1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 5, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + 
    s(id, wfpc.1, wfpc.2, wfpc.3,
      bs = "unc_pcre", xt = list("mfpc" = mfpca)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
sink("find_sim2.txt")
b_sim2 <- bamlss(f1, family = mjm_bamlss, data = d_indepri$data, 
                 timevar = "obstime", maxit = 1200, verbose_sampler = TRUE)
sink()

# Plot --------------------------------------------------------------------


set.seed(188)
ids <- sample(1:150, 5)
p <- ggplot(d_indepri$data %>% filter(id %in% ids) %>% 
         mutate(marker = recode_factor(marker, "m1" = "Marker 1",
                                       "m2" = "Marker 2")), 
       aes(x = obstime, y = y, color = id)) +
  theme_bw(base_size = 22) +
  geom_line(size = 1.2) +
  facet_grid(marker~.) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time")

pdf("../../../JointModel/Fig3.pdf", width = 6, height = 6)
p
dev.off()

times <- seq(0, 20, 1)
newdat <- d_indepri$data_short[rep(1:150, each = length(times)), ]
newdat$obstime <- rep(times, 150)
newdat <- newdat[newdat$survtime > newdat$obstime, ]
wfpcs <- eval_mfpc(mfpca, timepoints = newdat$obstime)
newdat <- rbind(newdat,
                newdat %>% mutate(marker = fct_recode(marker, "m2" = "m1")))
newdat <- cbind(newdat, wfpcs)
newdat$mu <- predict(b_sim, newdat, model = "mu")
newdat$fre <- predict(b_sim, newdat, model = "mu", term = "id", 
                      intercept = FALSE)


ggplot(data = newdat %>% filter(id %in% ids), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu), size = 1.2) +
  geom_line(data = d_indepri$data %>% filter(id %in% ids),
            aes(y = y), size = 0.5, linetype = "dashed") +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)

ggplot(newdat, aes(x = obstime, group = id)) + 
  geom_line(aes(y = wfpc.1), color = "red") + 
  geom_line(aes(y = wfpc.2), color = "blue") + 
  geom_line(aes(y = wfpc.3), color = "green") +
  facet_grid(~marker)

ggplot(newdat, aes(x = obstime, group = id, y = mu)) + 
  geom_line() +
  facet_grid(~marker)

newdat_id <- newdat %>% filter(id %in% ids) %>% droplevels() %>% split(.$id)
pars_id <- lapply(ids, function(i) (1+(i-1)*8):(1+(i-1)*8+7))
samp_means <- summary(b_sim$samples, quantiles = 0.5)$quantiles
newdat_id <- do.call(rbind, mapply(function(dat, pars) {
  ps <- samp_means[12:1211][pars]
  dat$fri_man <- c(as.matrix(dat[, 25:32]) %*% ps)
  dat
}, dat = newdat_id, pars = pars_id, SIMPLIFY = FALSE))

ggplot(data = newdat_id, aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fre), size = 1.2) +
  geom_line(aes(y = fri_man), size = 1.2, linetype = "dashed") +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)


ggplot(data = newdat %>% filter(id %in% ids), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu), size = 1.2) +
  geom_line(data = d_indepri$data %>% filter(id %in% ids),
            aes(y = y), size = 0.5, linetype = "dashed") +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(~marker)
