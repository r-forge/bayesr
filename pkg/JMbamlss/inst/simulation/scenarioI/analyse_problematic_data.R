
get_to_wd <- "../JMbamlss/"

source(paste0(get_to_wd, "inst/simulation/scenarioI/create_problematic_data.R"))
library(JMbayes2)
library(bamlss)
library(funData)
server_wd <- paste0("/run/user/1000/gvfs/smb-share:",
                    "server=clapton.wiwi.hu-berlin.de,",
                    "share=volkmana.hub/JMbamlss/simulation/")

# Convenience function
acc <- function(model){
  summ <- summary(model$samples[[1]][, grep("accepted", 
                                            colnames(model$samples[[1]]))])
  matrix(summ$statistics[, 1], ncol = 1, 
         dimnames = list(names(summ$statistics[, 1]), "Mean"))
}

# Data Plots --------------------------------------------------------------

ggplot(simdat,
       aes(x = year, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories",
          "Simulated Data Set 2")
plot(survfit(Surv(Time, event) ~ group, data = dat.id), 
     main = "Kaplan-Meier Estimates by 'Group'")



# Univariate Joint Model --------------------------------------------------

f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re") + s(year, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

# To recreate the problem, the first e.g. 350 observations are sufficient
set.seed(1)
b_uni6 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m6") %>% 
                   droplevels() %>%  
                   mutate(year = year + sqrt(.Machine$double.eps)) %>% 
                   as.data.frame(),
                 timevar = "year", idvar = "id", maxit = 1500, update.nu = TRUE,
                 verbose = TRUE) 
acc(b_uni6)





# # Data can also be censored at 10 with same problems
# simdat10 <- simdat %>%
#   mutate(event = ifelse(Time > 10, 0, event),
#          Time = ifelse(Time > 10, 10, Time)) %>%
#   filter(year <10) %>%
#   as.data.frame()
# ggplot(simdat10,
#        aes(x = year, y = y, colour = id)) +
#   facet_wrap(~marker, scales = "free_y") +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position = "none") +
#   ggtitle("Longitudinal Trajectories",
#           "Simulated Data Set 2")
# plot(survfit(Surv(Time, event) ~ group,
#              data = simdat10 %>% group_by(id) %>% slice_head(n = 1) %>% 
#                ungroup()),
#      main = "Kaplan-Meier Estimates by 'Group'")
# 
# set.seed(1)
# b_uni610 <- bamlss(f_uni, family = "jm", data = simdat10 %>%
#                      filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
#                    timevar = "year", idvar = "id", maxit = 1500,
#                    update.nu = TRUE, verbose = TRUE)
# summary(b_uni610$samples[[1]][, grep("accepted", 
#                                      colnames(b_uni610$samples[[1]]))])



# # Data can also be scaled by 10 with same problems
# simdat10 <- simdat %>%
#   mutate(y = y /10) %>%
#   as.data.frame()
# set.seed(1)
# b_uni610 <- bamlss(f_uni, family = "jm", data = simdat10 %>%
#                      filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
#                    timevar = "year", idvar = "id", maxit = 1500,
#                    update.nu = TRUE, verbose = TRUE)
# summary(b_uni610$samples[[1]][, grep("accepted", 
#                                      colnames(b_uni610$samples[[1]]))])

# Univariate JMbayes2 -----------------------------------------------------


# Get the quantile-based knots for comparability
kn <- mgcv::smoothCon(mgcv::s(Time, k = 20, bs = "ps"), 
                      data = dat)[[1]]$knots
CoxFit <- coxph(Surv(Time, event) ~ group, 
                data = dat.id)
lm6 <- lme(y6 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
set.seed(1)
jmb <- jm(CoxFit, list(lm6), time_var = "year",
          n_iter = 5500L, n_burnin = 500L, n_thin = 5L, 
          cores = 1, n_chains = 1, 
          GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
          save_random_effects = TRUE)

# Acceptance rates
str(jmb$acc_rates)
colMeans(jmb$acc_rates$b)



# Compare Samples ---------------------------------------------------------

# Bamlss Extract
rs_par <- grep("year,id", colnames(b_uni6$samples[[1]]))
rs_id <- rs_par[1:500]
rs_oth <- rs_par[501:504]


# True and estimated slopes
dat_tes <- data.frame(
  id = rep(seq_len(500), 3),
  rs = c(b[, 12], summary(b_uni6$samples[[1]][, rs_id])$statistics[, 1],
         jmb$statistics$Mean$b[, 2]),
  type = factor(rep(c(1:3), each = 500),
                labels = c("True", "bamlss", "JMbayes2"))
)
ggplot(dat_tes, aes(x = rs, colour = type)) +
  geom_density() +
  labs(colour = "", x = "Random Slope Parameters", y = "Density (Est + True)") +
  theme_bw() +
  stat_function(fun = dnorm, args = c(0, sqrt(0.197)), linetype = "dashed",
                color = "black") + 
  ggtitle("Empirical (and theoretical) Densities of RS Parameters",
          "Estimated based on Posterior Means/Simulated Parameters")


# Compare all draws with the variance parameters (BAMLSS)
all_draws <- data.frame(
  id = rep(seq_len(500), each = 1001),
  rs = c(b_uni6$samples[[1]][, rs_id])
)
rs_tau2 <- summary(b_uni6$samples[[1]][, rs_oth])$statistics[1, 1]
ggplot(all_draws, aes(x = rs)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  stat_function(fun = dnorm, args = c(0, sqrt(rs_tau2))) + 
  stat_function(fun = dnorm, args = c(0, sqrt(0.197)), col = "blue", 
                linetype = "dashed") + 
  theme_bw() +
  ggtitle("Distribution of All Random Slope Effect Samples", 
          "Univariate JM m6 (bamlss), 500 Ids, 1000 Samples, Truth (Blue)") +
  labs(x = "Random Slope Parameters", y = "Density")

# Compare all draws with the variance parameters
all_draws_jmb <- data.frame(
  id = rep(seq_len(500), 1000),
  rs = c(jmb$mcmc$b[[1]][, 2, ])
)
rs_tau2_jmb <- jmb$statistics$Mean$D[3]
ggplot(all_draws_jmb, aes(x = rs)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  stat_function(fun = dnorm, args = c(0, unname(sqrt(rs_tau2_jmb)))) + 
  stat_function(fun = dnorm, args = c(0, sqrt(0.197)), col = "blue", 
                linetype = "dashed") + 
  theme_bw() +
  ggtitle("Distribution of All Random Slope Effect Samples", 
          "Univariate JM m6 (JMbayes2), 500 Ids, 1000 Samples, Truth (Blue)") +
  labs(x = "Random Slope Parameters", y = "Density")


# Compare different Posteriori Hypers
rs_tau2_lo <- summary(b_uni6$samples[[1]][, rs_oth])$quantiles[1, 1]
rs_tau2_up <- summary(b_uni6$samples[[1]][, rs_oth])$quantiles[1, 5]
rs_tau2_lo_jmb <- jmb$statistics$CI_low$D[3]
rs_tau2_up_jmb <- jmb$statistics$CI_upp$D[3]


ggplot(data.frame(x = c(-3, 3)), aes(x = x)) +
  stat_function(fun = dnorm, args = c(0, sqrt(0.197)), col = "#F8766D", n = 10001) + 
  stat_function(fun = dnorm, args = c(0, sqrt(rs_tau2)), col = "#00BA38") + 
  stat_function(fun = dnorm, args = c(0, unname(sqrt(rs_tau2_jmb))), col = "#619CFF") + 
  stat_function(fun = dnorm, args = c(0, unname(sqrt(rs_tau2_lo_jmb))), col = "#619CFF",
                linetype = "dashed") + 
  stat_function(fun = dnorm, args = c(0, unname(sqrt(rs_tau2_up_jmb))), col = "#619CFF",
                linetype = "dashed") + 
  stat_function(fun = dnorm, args = c(0, sqrt(rs_tau2_lo)), col = "#00BA38",
                linetype = "dashed") + 
  stat_function(fun = dnorm, args = c(0, sqrt(rs_tau2_up)), col = "#00BA38",
                linetype = "dashed") + 
  theme_bw() +
  ggtitle("Normal Distributions with Posterior Variance of Random Slope Parameters", 
          "Simulation (Red), Posterior Mean/2.5%-Q/97.5%-Q: JMbayes2 (Blue), bamlss (Green)") +
  labs(x = "Random Slope Parameters", y = "Density") +
  annotate("text", x = -2, y = 0.875, parse = TRUE, label = "TrueTau2: 0.197") +
  annotate("text", x = -2, y = 0.625, parse = TRUE, label = paste0("JMBest:", 
                                                                  rs_tau2_jmb)) +
  annotate("text", x = 2, y = 0.75, parse = TRUE, label = paste0("BAMLSSest:", 
                                                                  rs_tau2))
rs_tau2




# Multivariate Joint Model ------------------------------------------------

# Basis functions on each dimension
seq1 <- seq(0, max(simdat$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs,
                                  scores = b)
mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })

# simdat <- simdat %>%
#   mutate(year = year + sqrt(.Machine$double.eps))
simdat_tru <- JMbamlss:::attach_wfpc(mfpca_tru, simdat,
                                     n = length(mfpca_tru$values),
                                     obstime = "year")
nfpc <- 12
f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1)
b_mul <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru, 
                timevar = "year", maxit = 1500, verbose  = TRUE,
                par_trace = TRUE)
# Crashes relaitvely late

set.seed(1)
b_mul <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru, 
                timevar = "year", maxit = 493, verbose  = TRUE, 
                sampler = FALSE, par_trace = TRUE)
saveRDS(b_mul, file = paste0(server_wd, "scen_mauff/mul/bamlss_mul_",
                             "nonuup.Rds"))
b_mul <- readRDS(paste0(server_wd, "scen_mauff/mul/bamlss_mul_",
                        "nonuup.Rds"))

# Multivariate JMbayes2 ---------------------------------------------------

lm1 <- lme(y1 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm2 <- lme(y2 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm3 <- lme(y3 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm4 <- lme(y4 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
lm5 <- lme(y5 ~ year, random = ~ year | id, 
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
mb <- jm(CoxFit, list(lm1, lm2, lm3, lm4, lm5, lm6), time_var = "year",
         n_iter = 1900L, n_burnin = 1000L, n_thin = 3L, 
         cores = 1, n_chains = 1, 
         GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
         save_random_effects = TRUE)