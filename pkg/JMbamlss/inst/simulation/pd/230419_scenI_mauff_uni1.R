
get_to_wd <- "../JMbamlss/"

source(paste0(get_to_wd, "inst/simulation/scenarioI/create_problematic_data.R"))
results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                     "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
library(JMbayes2)
library(bamlss)
library(funData)


# Univariate Joint Model --------------------------------------------------

f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  ti(year,bs="ps") + ti(id, bs = "re") + ti(year, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

d <- simdat %>%
                   filter(marker == "m1") %>% droplevels() %>%
                   as.data.frame()

d$year[abs(d$year) < 0.0001] <- 0.0001
d$Time[abs(d$Time) < 0.0001] <- 0.0001

set.seed(1)

b_uni1 <- bamlss(f_uni, family = "jm", data = d,
                 timevar = "year", idvar = "id", verbose = TRUE, update.nu = FALSE, maxit = 20, criterion = "BIC")

saveRDS(b_uni1, file = paste0(results_wd, "scen_mauff/uni/bamlss_uni1.Rds"))
b_uni1 <- readRDS(paste0(results_wd, "scen_mauff/uni/bamlss_uni1.Rds"))


# Univariate JMbayes2 -----------------------------------------------------

#Get the quantile-based knots for comparability
kn <- mgcv::smoothCon(mgcv::s(Time, k = 20, bs = "ps"),
                      data = dat)[[1]]$knots
CoxFit <- coxph(Surv(Time, event) ~ group,
                data = dat.id)
lm1 <- lme(y1 ~ year, random = ~ year | id,
           data = dat, na.action = na.omit,
           control = lmeControl(opt = "optim"))
set.seed(1)
jmb <- jm(CoxFit, list(lm1), time_var = "year",
          n_iter = 5500L, n_burnin = 500L, n_thin = 5L,
          cores = 1, n_chains = 1,
          GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
          save_random_effects = TRUE)
saveRDS(jmb, file = paste0(results_wd, "scen_mauff/uni/jmb_uni1.Rds"))
jmb <- readRDS(paste0(results_wd, "scen_mauff/uni/jmb_uni1.Rds"))


# Compare the Random Slopes -----------------------------------------------


# Bamlss Extract
rs_par <- grep("year,id", colnames(b_uni1$samples[[1]]))
rs_id <- rs_par[1:500]
rs_oth <- rs_par[501:504]


# True and estimated slopes
dat_tes <- data.frame(
  id = rep(seq_len(500), 3),
  rs = c(b[, 2], summary(b_uni1$samples[[1]][, rs_id])$statistics[, 1],
         jmb$statistics$Mean$b[, 2]),
  type = factor(rep(c(1:3), each = 500),
                labels = c("True", "bamlss", "JMbayes2"))
)
ggplot(dat_tes, aes(x = rs, colour = type)) +
  geom_density() +
  labs(colour = "", x = "Random Slope Parameters", y = "Density (Est + True)") +
  theme_bw() +
  stat_function(fun = dnorm, args = c(0, sqrt(0.006)), linetype = "dashed",
                color = "black") + 
  ggtitle("Empirical (and theoretical) Densities of RS Parameters",
          "Estimated based on Posterior Means/Simulated Parameters")


# Bamlss All Draws / tau2 -------------------------------------------------


# Compare all draws with the variance parameters
all_draws <- data.frame(
  id = rep(seq_len(500), each = 1000),
  rs = c(b_uni1$samples[[1]][, rs_id])
)
rs_tau2 <- summary(b_uni1$samples[[1]][, rs_oth])$statistics[1, 1]
ggplot(all_draws, aes(x = rs)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  stat_function(fun = dnorm, args = c(0, sqrt(rs_tau2))) + 
  stat_function(fun = dnorm, args = c(0, sqrt(0.006)), col = "blue", 
                linetype = "dashed") + 
  theme_bw() +
  ggtitle("Distribution of All Random Slope Effect Samples", 
          "Univariate JM m1 (bamlss), 500 Ids, 1000 Samples, Truth (Blue)") +
  labs(x = "Random Slope Parameters", y = "Density")



# JMbayes2 All Draws / tau2 -----------------------------------------------

# Compare all draws with the variance parameters
all_draws_jmb <- data.frame(
  id = rep(seq_len(500), 1000),
  rs = c(jmb$mcmc$b[[1]][, 2, ])
)
rs_tau2_jmb <- jmb$statistics$Mean$D[3]
ggplot(all_draws_jmb, aes(x = rs)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  stat_function(fun = dnorm, args = c(0, unname(sqrt(rs_tau2_jmb)))) + 
  stat_function(fun = dnorm, args = c(0, sqrt(0.006)), col = "blue", 
                linetype = "dashed") + 
  theme_bw() +
  ggtitle("Distribution of All Random Slope Effect Samples", 
          "Univariate JM m1 (JMbayes2), 500 Ids, 1000 Samples, Truth (Blue)") +
  labs(x = "Random Slope Parameters", y = "Density")



# Compare different Posteriori Hypers -------------------------------------


rs_tau2_lo <- summary(b_uni1$samples[[1]][, rs_oth])$quantiles[1, 1]
rs_tau2_up <- summary(b_uni1$samples[[1]][, rs_oth])$quantiles[1, 5]
rs_tau2_lo_jmb <- jmb$statistics$CI_low$D[3]
rs_tau2_up_jmb <- jmb$statistics$CI_upp$D[3]


ggplot(data.frame(x = c(-3, 3)), aes(x = x)) +
  stat_function(fun = dnorm, args = c(0, sqrt(0.006)), col = "#F8766D", n = 10001) + 
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
  labs(x = "Random Slope Parameters", y = "Density")


# Look at IWLS Steps ------------------------------------------------------

set.seed(1)
b_iwls <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m1") %>% droplevels() %>%
                   as.data.frame(),
                 timevar = "year", idvar = "id", optimizer = FALSE, 
                 start = parameters(b_uni1),
                 n.iter = 1200, burnin = 200, thin = 5, verbose = TRUE, 
                 save_iwls = TRUE)
IWLS_m1 <- IWLS 

mu_ri_compare <- sapply(IWLS_m1, function(it) {
  eq <- all.equal(it$mu[[1]]$mu0, it$mu[[1]]$mu1)
  if (is.logical(eq)) eq else FALSE
})

mu_rs_compare <- sapply(IWLS_m1, function(it) {
  eq <- all.equal(it$mu[[2]]$mu0, it$mu[[2]]$mu1, tolerance = 1e-14)
  if (is.logical(eq)) eq else FALSE
})

sigma_ri_compare <- sapply(IWLS_m1, function(it) {
  eq <- all.equal(it$mu[[1]]$xhess0 + it$mu[[1]]$phess0, 
                  it$mu[[1]]$xhess1 + it$mu[[1]]$phess1)
  if (is.logical(eq)) eq else FALSE
})
sigma_rs_compare <- sapply(IWLS_m1, function(it) {
  eq <- all.equal(it$mu[[2]]$xhess0 + it$mu[[2]]$phess0, 
                  it$mu[[2]]$xhess1 + it$mu[[2]]$phess1)
  if (is.logical(eq)) eq else FALSE
})
