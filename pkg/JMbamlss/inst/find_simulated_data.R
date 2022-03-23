
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

plot(d_indepri$data_short$survtime)
ggplot(d_indepri$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
