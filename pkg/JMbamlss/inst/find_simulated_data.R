
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
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 3),
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
sink("find_sim.txt")
b_sim <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                timevar = "obstime", verbose_sampler = TRUE)
sink()


set.seed(188)
ids <- sample(1:150, 5)
ggplot(d_indepri$data %>% filter(id %in% ids), 
       aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")

times <- seq(0, 20, 1)
newdat <- d_indepri$data_short[rep(1:150, each = length(times)), ]
newdat$obstime <- rep(times, 150)
newdat <- newdat[newdat$survtime > newdat$obstime, ]
wfpcs <- eval_mfpc(mfpca, timepoints = newdat$obstime)
newdat <- rbind(newdat,
                newdat %>% mutate(marker = fct_recode(marker, "m2" = "m1")))
newdat <- cbind(newdat, wfpcs)
newdat$mu <- predict(b_sim, newdat, model = "mu")

newdat <- left_join(newdat,
                    dat %>% mutate(y = b$y[[1]][, 3]) %>% 
                      select("id", "obstime", "y"),
                    by = c("id", "obstime"))
