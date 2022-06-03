library(survival)
library(bamlss)
library(MFPCA)
library(tidyverse)
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")


# Generate data with independent random intercepts
d_indepri <- simMultiJM(nsub = 150, times = seq(0, 1, by = 0.01), 
                        probmiss = 0.75, maxfac = 2,
                        nmark = 2, param_assoc = FALSE, M = 3, 
                        FPC_bases = NULL, 
                        FPC_evals = c(0.8, 0.5, 0.3),
                        mfpc_args = list(type = "split", eFunType = "PolyHigh",
                                         ignoreDeg = 1, eValType = "linear",
                                         eValScale = 1),
                        ncovar = 2,
                        lambda = function(t, x) {
                          1.6 * t^(0.6)
                        },
                        gamma = function(x) {
                          - 2 + 0.48*x[, 3]
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
                          log(0.2) + 0*t
                        }, 
                        tmax = NULL, seed = 1808, 
                        full = TRUE, file = NULL)
#save(d_indepri, file = "inst/objects/find_sim_data01.Rdata")
ggplot(d_indepri$data_short, aes(y = survtime, x = as.numeric(id), 
                                 color = factor(event))) +
  geom_point()
ggplot(d_indepri$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")

ggplot(d_indepri$data_full, aes(x = obstime, y = mu, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")


f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10, bs = "ps"),
  gamma ~ 1 + x3,
  mu ~ -1 + marker + obstime:marker + x3:marker + 
    s(id, fpc.1, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[1]])) + 
    s(id, fpc.2, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[2]])) + 
    s(id, fpc.3, bs = "unc_pcre", xt = list("mfpc" = mfpca_list[[3]])),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

mfpca <- create_true_MFPCA(M = 3, nmarker = 2, argvals = seq(0, 1, by = 0.01),
                           type = "split", eFunType = "PolyHigh",
                           ignoreDeg = 1, eValType = "linear",
                           eValScale = 1, evals = c(0.8, 0.5, 0.3))

mfpca_list <- lapply(1:3, function (i, m = mfpca) {
  list(functions = extractObs(m$functions, i),
       values = m$values[i])
})


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
compile_alex("~/Dokumente/Arbeit/Greven/JM/JMbamlss/src")

set.seed(1808)
sink("find_sim01.txt")
b_sim1 <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                timevar = "obstime", maxit = 1200, verbose_sampler = TRUE)
sink()
save(b_sim1, file = "inst/objects/find_sim01.Rdata")


# # Why is sigma not accepted? ----------------------------------------------
# 
# 
# b <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
#             timevar = "obstime", optimizer = FALSE, start = parameters(b_sim1))
# mfpca <- create_true_MFPCA(M = 3, nmarker = 2, argvals = seq(0, 25, by = 0.25),
#                            type = "split", eFunType = "PolyHigh",
#                            ignoreDeg = 1, eValType = "linear",
#                            eValScale = 1, evals = c(80:78))
# f1 <- list(
#   Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 5, bs = "ps"),
#   gamma ~ 1 + x3,
#   mu ~ -1 + marker + obstime:marker + x3:marker + 
#     s(id, wfpc.1, wfpc.2, wfpc.3,
#       bs = "unc_pcre", xt = list("mfpc" = mfpca)),
#   sigma ~ -1 + marker,
#   alpha ~ -1 + marker
# )
# sink("find_sim2.txt")
# b_sim2 <- bamlss(f1, family = mjm_bamlss, data = d_indepri$data, 
#                  timevar = "obstime", maxit = 1200, verbose_sampler = TRUE)
# sink()

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
newdat <- d_indepri$data_full
# fpcs <- eval_mfpc(mfpca, timepoints = newdat$obstime, marker = newdat$marker)
# newdat <- rbind(newdat,
#                 newdat %>% mutate(marker = fct_recode(marker, "m2" = "m1")))
# newdat <- cbind(newdat, wfpcs)
newdat$mu_pred <- predict(b_sim1, newdat, model = "mu")
newdat$fre <- predict(b_sim1, newdat, model = "mu", term = "id", 
                      intercept = FALSE)


ggplot(data = newdat %>% filter(id %in% ids), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu), size = 1.2) +
  geom_line(aes(y = mu_pred), size = 1.2, linetype = "dotted") +
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
pars_id <- lapply(sort(ids), function(i) c(18+i, 172+i, 326+i))
samp_means <- summary(b_sim1$samples, quantiles = 0.5)$quantiles
newdat_id <- do.call(rbind, mapply(function(dat, pars) {
  ps <- samp_means[pars]
  dat$fri_man <- c(as.matrix(dat[, 20:22]) %*% ps)
  dat
}, dat = newdat_id, pars = pars_id, SIMPLIFY = FALSE))

ggplot(data = newdat_id, aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fre), linetype = "dotted") +
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
