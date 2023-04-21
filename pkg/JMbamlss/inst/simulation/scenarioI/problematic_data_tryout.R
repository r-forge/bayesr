
acc <- function(model){
  summ <- summary(model$samples[[1]][, grep("accepted", 
                                             colnames(model$samples[[1]]))])
  matrix(summ$statistics[, 1], ncol = 1, 
         dimnames = list(names(summ$statistics[, 1]), "Mean"))
}

# Analyze The Posterior Mode ----------------------------------------------


# Traceplots
lambda <- sapply(b_uni6$model.stats$optimizer$par_trace, 
                 function(x) x$lambda$s[[1]])
gamma <- sapply(b_uni6$model.stats$optimizer$par_trace, 
                function(x) x$gamma$p)
alpha <- sapply(b_uni6$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
mus1 <- sapply(b_uni6$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[1]])
mus2 <- sapply(b_uni6$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[2]])
mup <- sapply(b_uni6$model.stats$optimizer$par_trace, 
              function(x) x$mu$p)
sigma <- sapply(b_uni6$model.stats$optimizer$par_trace, 
                function(x) x$sigma$p)

ggplot(data = data.frame(t(lambda)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Parameters")

ggplot(data = data.frame(t(gamma)) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -5.8, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Gamma Parameters")

ggplot(data = data.frame(alpha) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")

ggplot(data = data.frame(t(mus1)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Intercept Parameters")

ggplot(data = data.frame(t(mus2)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Slope Parameters")

ggplot(data = data.frame(t(mup)) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -0.3, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 1.71, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Parameters")

ggplot(data = data.frame(sigma) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = log(0.57), linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Sigma Parameters")


# Estimated Baseline Hazards
lambga_dat <- data.frame("Time" = seq(0, 30, length.out = 100))
x_lambga <- mgcv::smoothCon(mgcv::s(Time, k = 20, bs = "ps"), data = lambga_dat,
                            absorb.cons = TRUE)[[1]]$X
lambga_dat$fit <- x_lambga %*% b_uni6$parameters$lambda$s$`s(Time)`[1:19] + 
  b_uni6$parameters$gamma$p[["(Intercept)"]]
lambga_dat$tru <- log(1.65) + 0.65*log(lambga_dat$Time) - 5.8
ggplot(lambga_dat, aes(x = Time)) +
  geom_line(aes(y = fit)) +
  geom_line(aes(y = tru), linetype = "dotted") +
  theme_bw() +
  labs(y = "Lambda + Gamma Intercept") +
  ggtitle("Baseline Hazard Parameters")


# Estimated Random Effects
b_ri <- b[1:350, 11]
b_ri_dif <- sapply(seq_len(350), function (i) b_ri[i] - mus1[i,])
ggplot(data = data.frame(b_ri_dif) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "b - b_hat") +
  ggtitle("Mu Random Intercept Bias")

b_ri_bias_dif <- data.frame(b_ri_dif) %>%
  mutate(it = 1:50) %>%
  pivot_longer(cols = -it) %>% 
  arrange(name, it) %>%
  group_by(name) %>%
  mutate(DIFF = value - lag(value)) %>%
  as.data.frame() %>% 
  na.omit()
ggplot(b_ri_bias_dif, aes(x = it, y = DIFF, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "(b - b_hat)_n - (b - b_hat)_(n-1)") +
  ggtitle("Mu Random Intercept Bias Changes")



b_rs <- b[1:350, 12]
b_rs_dif <- sapply(seq_len(350), function (i) b_rs[i] - mus2[i,])
ggplot(data = data.frame(b_rs_dif) %>%
         mutate(it = 1:50) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "b - b_hat") +
  ggtitle("Mu Random Slope Bias")

b_rs_bias_dif <- data.frame(b_rs_dif) %>%
  mutate(it = 1:50) %>%
  pivot_longer(cols = -it) %>% 
  arrange(name, it) %>%
  group_by(name) %>%
  mutate(DIFF = value - lag(value)) %>%
  as.data.frame() %>% 
  na.omit()
ggplot(b_rs_bias_dif, aes(x = it, y = DIFF, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "(b - b_hat)_n - (b - b_hat)_(n-1)") +
  ggtitle("Mu Random Slope Bias Changes")



# Change the Sampling Step-Length -----------------------------------------

set.seed(1)
b_0.1 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                timevar = "year", idvar = "id", optimizer = FALSE,
                start = parameters(b_uni6), verbose = TRUE, nu_sam = 0.1) 
summary(b_0.1$samples[[1]][, grep("accepted", colnames(b_0.1$samples[[1]]))])


set.seed(1)
b_2 <- bamlss(f_uni, family = "jm", data = simdat %>%
                  filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                timevar = "year", idvar = "id", optimizer = FALSE,
                start = parameters(b_uni6), verbose = TRUE, nu_sam = 2) 
summary(b_2$samples[[1]][, grep("accepted", colnames(b_2$samples[[1]]))])



# Without Burnin ----------------------------------------------------------

set.seed(1)
b_nobu2 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                 timevar = "year", idvar = "id", optimizer = FALSE,
                 start = parameters(b_uni6), verbose = TRUE, n.iter = 500,
                 burnin = 0, thin = 1) 
acc(b_nobu)

# Change the bamlss Code and shift the tau update to beginning of 
# propose_mu_matrix
# set.seed(1)
# b_nobu2 <- bamlss(f_uni, family = "jm", data = simdat %>%
#                    filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
#                  timevar = "year", idvar = "id", optimizer = FALSE,
#                  start = parameters(b_uni6), verbose = TRUE, n.iter = 500,
#                  burnin = 0, thin = 1) 
# acc(b_nobu2)

# Change the prior assumption to half Cauchy
f_hc <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps", 
                                       xt = list("prior" = "hc")),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re", xt = list("prior" = "hc")) + 
    s(year, id, bs = "re", xt = list("prior" = "hc")),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
set.seed(1)
b_nbhc <- bamlss(f_hc, family = "jm", data = simdat %>%
                   filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                 timevar = "year", idvar = "id", optimizer = FALSE,
                 start = parameters(b_uni6), verbose = TRUE, n.iter = 500,
                 burnin = 0, thin = 1) 
acc(b_nbhc)

# Change the prior assumption to half Normal
f_hn <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps",
                                       xt = list("prior" = "hn")),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re", xt = list("prior" = "hn")) + 
    s(year, id, bs = "re", xt = list("prior" = "hn")),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
set.seed(1)
b_nbhn <- bamlss(f_hn, family = "jm", data = simdat %>%
                   filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                 timevar = "year", idvar = "id", optimizer = FALSE,
                 start = parameters(b_uni6), verbose = TRUE, n.iter = 500,
                 burnin = 0, thin = 1) 
acc(b_nbhn)


# Save IWLS Info ----------------------------------------------------------

set.seed(1)
b_iwls6 <- bamlss(f_uni, family = "jm", data = simdat %>%
                    filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                  timevar = "year", idvar = "id", optimizer = FALSE,
                  start = parameters(b_uni6), verbose = TRUE, save_iwls = TRUE)
IWLS_m6 <- IWLS
acc(b_iwls6)
IWLS_m6

mu_ri_compare6 <- sapply(IWLS_m6, function(it) {
  eq <- all.equal(it$mu[[1]]$mu0, it$mu[[1]]$mu1)
  if (is.logical(eq)) eq else FALSE
})

mu_rs_compare6 <- sapply(IWLS_m6, function(it) {
  eq <- all.equal(it$mu[[2]]$mu0, it$mu[[2]]$mu1, tolerance = 1e-14)
  if (is.logical(eq)) eq else FALSE
})

sigma_ri_compare6 <- sapply(IWLS_m6, function(it) {
  eq <- all.equal(it$mu[[1]]$xhess0 + it$mu[[1]]$phess0, 
                  it$mu[[1]]$xhess1 + it$mu[[1]]$phess1)
  if (is.logical(eq)) eq else FALSE
})
sigma_rs_compare6 <- sapply(IWLS_m6, function(it) {
  eq <- all.equal(it$mu[[2]]$xhess0 + it$mu[[2]]$phess0, 
                  it$mu[[2]]$xhess1 + it$mu[[2]]$phess1)
  if (is.logical(eq)) eq else FALSE
})
