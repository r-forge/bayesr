
get_to_wd <- "../JMbamlss/"

source(paste0(get_to_wd, "inst/simulation/scenarioI/create_problematic_data.R"))
library(JMbayes2)
library(bamlss)
library(funData)

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
                   filter(marker == "m6", id %in% 1:350) %>% as.data.frame(),
                 timevar = "year", idvar = "id", maxit = 1500, update.nu = TRUE,
                 verbose = TRUE, par_trace = TRUE) 
summ <- summary(b_uni6$samples[[1]][, grep("accepted", 
                                           colnames(b_uni6$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))


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
jmb <- jm(CoxFit, list(lm6), time_var = "year",
          n_iter = 1900L, n_burnin = 1000L, n_thin = 3L, 
          cores = 1, n_chains = 1, 
          GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
          save_random_effects = TRUE)

# Acceptance rates
str(jmb$acc_rates)
colMeans(jmb$acc_rates$b)


# Multivariate Joint Model ------------------------------------------------

# Basis functions on each dimension
seq1 <- seq(0, max(simdat$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs)
mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })

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
b <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru, 
            timevar = "year", maxit = 1500, update_nu = TRUE, verbose  = TRUE)
# Crashes relaitvely late

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