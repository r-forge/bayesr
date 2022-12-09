library(tidyverse)
library(survival)

# Number of individuals and other quantities
n <- 300
argvals <- seq(0, 1, by = 0.01)
x <- seq(0, 1, by = 0.1)

# Random covariance matrix
# Set the eigenvalues but the eigenvectors are random
set.seed(1105)
p <- 6
P <- qr.Q(qr(matrix(rnorm(p^2), ncol = p)))
evals <- c(4, 3, 2, 1, 0.5, 0.2)
cov <- crossprod(P, P*(evals))


# Find spline functions
# Marker1
m1sp1 <- splinefun(x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
m1sp2 <- splinefun(x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
m1sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4))
# Marker2
m2sp1 <- splinefun(x, c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
m2sp2 <- splinefun(x, c(0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
m2sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3))

m1 <- funData(argvals = argvals,
              X = matrix(c(m1sp1(argvals), m1sp2(argvals), m1sp3(argvals)), 
                         nrow = 3, byrow = TRUE))
m2 <- funData(argvals = argvals,
              X = matrix(c(m2sp1(argvals), m2sp2(argvals), m2sp3(argvals)), 
                         nrow = 3, byrow = TRUE))
plot(m1)
title("Basis Functions on Marker1")
plot(m2)
title("Basis Functions on Marker2")

# True multivariate covariance structure
m <- MFPCA_cov(cov = cov, basis_funs = list(m1, m2))

# True eigenvalues
mevals <- matrix(c(m$values, cumsum(m$values)/sum(m$values)), nrow = 2,
                 dimnames = list(rownames = c("Mevals", "Cumsum")),
                 byrow = TRUE)
uni_norms <- lapply(m$functions, norm)
uni_ratio <- t(sapply(uni_norms, function (n) {
  cumsum(m$values*n)/sum(m$values*n)
}))
xtable::xtable(rbind(mevals, uni_ratio))
ggplot(multifammPaper::funData2DataFrame(m$functions),
       aes(x = t, y = y, colour = factor(obs))) +
  geom_line() +
  facet_grid(dim~., labeller = as_labeller(c("1" = "Marker1", "2" ="Marker2"))) +
  theme_bw() +
  ggtitle("Multivariate Eigenfunctions") +
  labs(colour = "MFPC")


# Random draws from the multivariate Gaussian Random Process
n_dr <- 20
sc <- mvtnorm::rmvnorm(n = n_dr, sigma = diag(m$values))
mfData <- multiFunData(lapply(m$functions, function(x) {
  funData(argvals = argvals, X = sc %*% x@X)
}))
ggplot(multifammPaper::funData2DataFrame(mfData),
       aes(x = t, y = y, colour = factor(obs))) +
  geom_line() +
  facet_grid(dim~., labeller = as_labeller(c("1" = "Marker1", "2" ="Marker2"))) +
  theme_bw() +
  ggtitle("Random Observations") +
  theme(legend.position = "none")



# Simulate the data -------------------------------------------------------


d_rirs <- simMultiJM(nsub = n, times = seq(0, 1, by = 0.01), 
                     max_obs = 15, probmiss = 0.75, maxfac = 3,
                     nmark = 2, long_assoc = "FPC", M = 6, 
                     FPC_bases = m$functions, FPC_evals = m$values, 
                     #mfpc_args = NULL, re_cov_mat = NULL,
                     ncovar = 2,
                     lambda = function(t, x) {
                       1.65 * t^(0.65)
                     },
                     gamma = function(x) {
                       -3 + 0.3*x[, 3]
                     },
                     alpha = list(function(t, x) {
                       1.1 + 0*t
                     }, function(t, x) {
                       1.1 + 0*t
                     }),
                     mu = list(function(t, x, r){
                       0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                     }, function(t, x, r){
                       0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                     }),
                     sigma = function(t, x) {
                       log(0.06) + 0*t
                     }, 
                     tmax = NULL, seed = 1008, 
                     full = TRUE, file = NULL)


ggplot(data = d_rirs$data_short %>% slice_head(n = n),
       aes(x = id, y = survtime, colour = as.factor(event))) +
  geom_point()
ggplot(d_rirs$data, aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
ggplot(d_rirs$data_hypo %>% filter(event == 1), aes(x = obstime, y = mu, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none")
summary(c(table(d_rirs$data$id)))
table(d_rirs$data_short[seq_len(n), "event"])



# Large Sample ------------------------------------------------------------

n_high <- 30000
d_rirs <- simMultiJM(nsub = n_high, times = seq(0, 1, by = 0.01), 
                     max_obs = 15, probmiss = 0.75, maxfac = 3,
                     nmark = 2, long_assoc = "FPC", M = 6, 
                     FPC_bases = m$functions, FPC_evals = m$values, 
                     #mfpc_args = NULL, re_cov_mat = NULL,
                     ncovar = 2,
                     lambda = function(t, x) {
                       1.65 * t^(0.65)
                     },
                     gamma = function(x) {
                       -3 + 0.3*x[, 3]
                     },
                     alpha = list(function(t, x) {
                       1.1 + 0*t
                     }, function(t, x) {
                       1.1 + 0*t
                     }),
                     mu = list(function(t, x, r){
                       0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                     }, function(t, x, r){
                       0 + 1*t + 0.3*x[, 3] + 0.3*t*x[, 3]
                     }),
                     sigma = function(t, x) {
                       log(0.06) + 0*t
                     }, 
                     tmax = NULL, seed = 1008, 
                     full = TRUE, file = NULL)
table(d_rirs$data_short[seq_len(n_high), "event"]) /n_high
ggplot(d_rirs$data_short[seq_len(n_high), ], 
       aes(x = survtime, colour = as.factor(event) )) +
  geom_density()+
  ggtitle("Density of Follow-Up Time (Censoring or Event)") +
  theme_bw() +
  labs(colour = "Event", y = "Density", x = "Follow-Up Time")
fit1 <- survfit(Surv(survtime, event) ~ x3, 
                data = d_rirs$data_short[seq_len(n_high), ])
plot(fit1)
title("Kaplan Meier Plot based on X3")
ggplot(d_rirs$data %>% filter(id %in% seq_len(100)),
       aes(x = obstime, y = y, color = id)) +
  geom_line() +
  facet_grid(~marker) +
  theme(legend.position = "none") +
  ggtitle("Longitudinal Trajectories for 100 Random Draws")
ggplot(d_rirs$data_hypo %>% filter(id %in% seq_len(100)),
       aes(x = obstime, y = mu, color = id)) +
  geom_line() +
  facet_grid(event~marker) +
  theme(legend.position = "none")+
  ggtitle("Hypothetical Longitudinal Trajectories for 100 Random Draws By Event")
