
# Set Up R Session --------------------------------------------------------

library(MFPCA)
library(tidyverse)
library(Matrix)
library(mvtnorm)
library(survival)
library(bamlss)


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
source("R/survint.R")
source("R/compile.R")
source("R/MJM_predict.R")

# Compile the C function
compile_alex()

d_indepri <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                        probmiss = 0.75, maxfac = 1.5,
                        nmark = 2, param_assoc = TRUE, M = NULL, 
                        FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                        re_cov_mat = matrix(c(0.68, 0, 0, 0.68), ncol = 2), 
                        ncovar = 2,
                        lambda = function(t, x) {
                          1.37 * t^(0.37)
                        },
                        gamma = function(x) {
                          - 5.8 + 0.48*x[, 3]
                        },
                        alpha = list(function(t, x) {
                          0.64 + 0*t
                        }, function(t, x) {
                          -0.64 + 0*t
                        }),
                        mu = list(function(t, x, r){
                          2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] + 
                            r[, 1]
                        }, function(t, x, r){
                          2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] + 
                            r[, 2]
                        }),
                        sigma = function(t, x) {
                          log(0.6) + 0*t
                        }, 
                        tmax = NULL, seed = 1808, 
                        full = TRUE, file = NULL)

load("inst/objects/indepri_models.Rdata")
summary(rowMeans(b_indepri_byre$samples[[1]][, c(19:168)]))
summary(rowMeans(b_indepri_byre$samples[[1]][, c(173:322)]))
summary(rowMeans(b_indepri_byre$samples[[1]][, c(19:322)]))


load("inst/objects/indeprirs_models.Rdata")
summary(rowMeans(b_indeprirs_re$samples[[1]][, c(19:318)]))
summary(rowMeans(b_indeprirs_re$samples[[1]][, c(323:622)]))
mean(rowMeans(b_indeprirs_re$samples[[1]][, c(19:168)]))
mean(rowMeans(b_indeprirs_re$samples[[1]][, c(169:318)]))
summary(b_indeprirs_re, model = "mu")
group0 <- (19:168)[d_indepri$data_short$x3[1:150] == 0]
group1 <- (19:168)[d_indepri$data_short$x3[1:150] == 1]
mean(rowMeans(b_indeprirs_re$samples[[1]][, group0]))
mean(rowMeans(b_indeprirs_re$samples[[1]][, group1]))
2.189442 + 0.003555687
2.189442 + 0.003555687 +
colnames(b_indeprirs_re$samples[[1]])[c(19, 168, 169, 318)]
colnames(b_indeprirs_re$samples[[1]])[c(323, 472, 473, 622)]


################################################################

b1 <- gam(score ~ Machine + s(Worker, bs = "re") + 
            s(Worker, bs = "re", by = Machine), data = Machines, method = "REML")
sum(b1$coefficients[4:9])
sum(b1$coefficients[10:15])
sum(b1$coefficients[16:21])
sum(b1$coefficients[22:27])


################################################################
set.seed(1808)
f <- list(
  score ~ Machine + s(Worker, bs = "re") + 
    s(Worker, bs = "re", by = Machine)
)
bamlss1 <- bamlss(f, data = Machines)
summary(bamlss1)
summary(rowSums(bamlss1$samples[[1]][, 1:6]))
summary(rowSums(bamlss1$samples[[1]][, 11:16]))
summary(rowSums(bamlss1$samples[[1]][, 21:26]))

################################################################

set.seed(1808)
n <- 100
b <- mvtnorm::rmvnorm(n = n, mean = c(0,0), sigma = diag(c(0.68, 0.28)))
t <- seq(0, 25, by = 0.25)
x <- sample(c(0, 1), size = n, replace = TRUE)
dat <- data.frame(id = factor(rep(seq_len(n), each = length(t))),
                  t = rep(t, n),
                  x = rep(x, each = length(t)))
dat$mu <- c(sapply(seq_len(n), function(i) {
  2.13 + 0.24*t - 0.25*x[i] - 0.05*t*x[i] + b[i, 1] + b[i, 2]*t
}))
n_i <- sample(15:40, size = n, replace = TRUE)
t_i <- lapply(n_i, function (i){
  sort(sample(x = seq_along(t), size = i, replace = FALSE))
})
dat <- do.call(rbind, mapply(function (d, ts) {
  d[ts, , drop = FALSE]
}, d = split(dat, dat$id), ts = t_i, SIMPLIFY = FALSE))
dat$y <- rnorm(n = nrow(dat), mean = dat$mu, sd = 0.6)

m1 <- bam(y ~ x + t + x:t + s(id, bs = "re") + s(id, t, bs = "re"), 
          data = dat, method = "fREML")
summary(m1)
sum(m1$coefficients[5:(n+4)])
sum(m1$coefficients[(n+5):(2*n+4)])

m2 <- bamlss(list(y ~ x + t + x:t + s(id, bs = "re") + s(id, t, bs = "re"),
                  sigma ~ 1),
             family = "gaussian", data = dat)
summary(m2)
newdat <- data.frame(id = factor(seq_len(n)),
                     t = rep(1, n))
ri <- predict(m2, model = "mu", term = "s(id)", newdata = newdat, 
              intercept = FALSE)
rs <- predict(m2, model = "mu", term = "s(id,t)", newdata = newdat,
              intercept = FALSE)
b[1,]
colSums(b)
summary(rowMeans(m2$samples[[1]][, 1:100]))
summary(rowSums(m2$samples[[1]][, 1:100]))

##################################################################

d_ri_onem <- simMultiJM(nsub = 150, times = seq(0, 25, by = 0.25), 
                        probmiss = 0.75, maxfac = 1.5,
                        nmark = 1, param_assoc = TRUE, M = NULL, 
                        FPC_bases = NULL, FPC_evals = NULL, mfpc_args = NULL,
                        re_cov_mat = matrix(c(0.68), ncol = 1), 
                        ncovar = 2,
                        lambda = function(t, x) {
                          1.37 * t^(0.37)
                        },
                        gamma = function(x) {
                          - 5.8 + 0.48*x[, 3]
                        },
                        alpha = list(function(t, x) {
                          0.64 + 0*t
                        }),
                        mu = list(function(t, x, r){
                          2.13 + 0.24*t - 0.25*x[, 3] - 0.05*t*x[, 3] +  r
                        }),
                        sigma = function(t, x) {
                          log(0.6) + 0*t
                        }, 
                        tmax = NULL, seed = 1808, 
                        full = TRUE, file = NULL)
ggplot(d_ri_onem$data, aes(x = obstime, color = id)) +
  geom_point(aes(y = y)) +
  geom_line(aes(y = mu)) +
  theme(legend.position = "none")

f_ri_onem <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 10),
  gamma ~ 1 + x3,
  mu ~ obstime + x3 + obstime:x3 + s(id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1808)
sink("ri_onem.txt")
b_ri_onem <- bamlss(f_ri_onem, family = "jm", data = d_ri_onem$data,
                    timevar = "obstime", idvar = "id", maxit = 1000)
sink()