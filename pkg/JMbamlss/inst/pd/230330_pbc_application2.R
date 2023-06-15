# PBC Data Analysis -------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


# Always
library(survival)
library(JMbayes2)
library(JMbamlss)
library(tidyverse)


# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
# Logarithmic scale for all three markers
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years + sqrt(.Machine$double.eps),
         obstime = year + sqrt(.Machine$double.eps), marker = factor(marker),
         logy = log(y)) %>%
  select(id, survtime, event, sex, drug, age, marker, obstime, y, logy) %>%
  arrange(marker, id, obstime) %>%
  na.omit() %>%
  as.data.frame()

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(dplyr::n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis)

ggplot(p_long, aes(x = obstime, y = y, group = id)) +
  geom_line() +
  facet_wrap(~marker, scales = "free") +
  theme_bw() +
  ggtitle("PBC Data")
ggplot(p_long, aes(x = obstime, y = logy, group = id)) +
  geom_line() +
  facet_wrap(~marker, scales = "free") +
  theme_bw() +
  ggtitle("PBC Data on Log Scale")


# Univariate Models -------------------------------------------------------

# Random intercept + random slope
f_uni <- list(
  Surv2(survtime, event, obs = y) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  mu ~ sex + s(obstime, xt = list("scale" = FALSE)) + 
    s(id, bs = "re", xt = list("scale" = FALSE)) + 
    s(obstime, id, bs = "re", xt = list("scale" = FALSE)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1559)
b_bil <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE,
                update.nu = TRUE)
b_cho <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "serChol"),
                timevar = "obstime", idvar = "id", maxit = 5000,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_sgo <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "SGOT"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

saveRDS(b_bil, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_uni_bil.Rds"))
saveRDS(b_cho, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_uni_cho.Rds"))
saveRDS(b_sgo, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_uni_sgo.Rds"))

# Random intercept + random slope log scale
f_uni_log <- list(
  Surv2(survtime, event, obs = logy) ~ -1 +
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  mu ~ sex + s(obstime, xt = list("scale" = FALSE)) + 
    s(id, bs = "re", xt = list("scale" = FALSE)) + 
    s(obstime, id, bs = "re", xt = list("scale" = FALSE)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1559)
b_bil <- bamlss(f_uni_log, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 5000)#,
                # n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE,
                # update.nu = TRUE)
b_cho <- bamlss(f_uni_log, family = "jm",
                data = p_long %>% filter(marker == "serChol"),
                timevar = "obstime", idvar = "id", maxit = 5000,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_sgo <- bamlss(f_uni_log, family = "jm",
                data = p_long %>% filter(marker == "SGOT"),
                timevar = "obstime", idvar = "id", maxit = 5000,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

saveRDS(b_bil, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_unilog_bil.Rds"))
saveRDS(b_cho, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_unilog_cho.Rds"))
saveRDS(b_sgo, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_unilog_sgo.Rds"))

summ <- summary(b_bil$samples[[1]][, grep("\\.accepted", 
                                          colnames(b_bil$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))

post_tau <- summary(b_cho$samples[[1]][, 
  grep("obstime,id).tau", colnames(b_cho$samples[[1]]))])$statistics[1]
post_re <- data.frame("rs" = summary(b_cho$samples[[1]][,
  grep("obstime,id)\\.obstime", colnames(b_cho$samples[[1]]))])$statistics[, 1])
ggplot(post_re, aes(x = rs)) +
  geom_density(col = "blue") +
  stat_function(fun = dnorm, args = c(0, unname(sqrt(post_tau))), col = "red") + 
  theme_bw()


ri_id <- grep("\\(id\\)", colnames(b_bil$samples[[1]]))[1:304]
rs_id <- grep("\\(obstime,id\\)", colnames(b_bil$samples[[1]]))[1:304]
all_draws <- data.frame(
  id = rep(rep(seq_len(304), each = 1001), times = 2),
  re = c(b_bil$samples[[1]][, ri_id], b_bil$samples[[1]][, rs_id]),
  type = factor(rep(c("Random Intercept", "Random Slope"), each = 1001*304))
)
ri_tau2 <- unname(summary(b_bil$samples[[1]][, ri_id[304] + 1])$statistics[1])
rs_tau2 <- unname(summary(b_bil$samples[[1]][, rs_id[304] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-2, 2, by = 0.01), seq(-2, 2, by = 0.01)),
  y = c(dnorm(seq(-2, 2, by = 0.01), 0, sqrt(ri_tau2)),
        dnorm(seq(-2, 2, by = 0.01), 0, sqrt(rs_tau2))),
  type = factor(c(rep("Random Intercept", length(seq(-2, 2, by = 0.01))),
                  rep("Random Slope", length(seq(-2, 2, by = 0.01)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("Distribution of All Random Slope Effect Samples", 
          paste0("Univariate JM pbc (bamlss), 304 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Effect Parameters", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")



# Change the linear predictors so that it corresponds to the code on the 
# homepage
f_uni_log <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + s(survtime, k = 20),
  gamma ~ 1 + s(age, k = 20) + drug + sex,
  mu ~ ti(obstime, k = 20) + ti(id, bs = "re") + ti(obstime, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
set.seed(1559)
b_bil <- bamlss(f_uni_log, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                verbose = TRUE, update.nu = TRUE)
saveRDS(b_bil, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_uniloghome_bil.Rds"))
summ <- summary(b_bil$samples[[1]][, grep("\\.accepted", 
                                          colnames(b_bil$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))




# Random intercept + random slope log scale + flexible alpha
f_uni_log_alpha <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  mu ~ sex + s(obstime) + s(id, bs = "re") + s(obstime, id, bs = "re"),
  sigma ~ 1,
  alpha ~ s(survtime, k = 20),
  dalpha ~ -1
)
set.seed(1559)
b_bil_alpha <- bamlss(f_uni_log_alpha, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 5000,
                n.iter = 500, burnin = 0, thin = 1, verbose = TRUE)
saveRDS(b_bil_alpha, file = paste0("~/Documents/joint_models/JointModel/",
                                   "PBC_analysis/pbc_unilog_bil_alpha.Rds"))
# Converges at 4295
set.seed(1559)
b_bil_alpha1 <- bamlss(f_uni_log_alpha, family = "jm",
                      data = p_long %>% filter(marker == "serBilir"),
                      timevar = "obstime", idvar = "id", optimizer = FALSE,
                      start = parameters(b_bil_alpha),
                      n.iter = 1500, burnin = 500, thin = 1, verbose = TRUE)
# 0 acceptance probabilities for random intercept


# Univariate Functional Random Intercept ----------------------------------


f_fri <- list(
  Surv2(survtime, event, obs = logy) ~ s(survtime, k = 20, bs = "ps", 
                                         xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  mu ~ sex + ti(id, bs = "re", xt = list("scale" = FALSE)) + 
    ti(obstime, xt = list("scale" = FALSE)) + 
    ti(id, obstime, bs = c("re", "cr"), k = c(nlevels(pbc2$id), 8),
       xt = list("scale" = FALSE)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1559)
b_bil <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 1500)#,
                # n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_cho <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "serChol") %>%
                  droplevels() %>% as.data.frame(),
                timevar = "obstime", idvar = "id", maxit = 5000)#,
                #n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_sgo <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "SGOT"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

# Add random noise
set.seed(1559)
b_bil2 <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "serBilir") %>%
                  mutate(logy = logy + rnorm(nrow(.))) %>%
                  droplevels() %>% as.data.frame(),
                timevar = "obstime", idvar = "id", maxit = 1500)


# Linear predictors as in website
f_fri <- list(
  Surv2(survtime, event, obs = logy) ~ s(survtime, k = 20),
  gamma ~ s(age, k = 20) + drug + sex,
  mu ~ ti(id, bs="re") + 
    ti(obstime, k=20) + 
    ti(id, obstime, bs = c("re", "cr"), k = c(nlevels(pbc2$id), 8)),
  sigma ~ 1,
  alpha ~ s(survtime, k = 20),
  dalpha ~ -1
)

set.seed(1505)
b_cho <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "serChol"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                verbose = TRUE)
set.seed(1505)
b_sgo <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "SGOT"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                verbose = TRUE)
saveRDS(b_cho, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_unilog_cho_fri.Rds"))
saveRDS(b_sgo, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_unilog_sgo_fri.Rds"))



# Univariate Models as calculated in bamlss-Homepage ----------------------

data("pbc2", package = "JMbayes")
## Set up the model formula including
## functional random intercepts using ti().
f <- list(
  Surv2(years, status2, obs = log(serBilir)) ~ s(years,k=20),
  gamma ~ s(age,k=20) + drug + sex,
  mu ~ ti(id,bs="re") + 
    ti(year,k=20) + 
    ti(id,year,bs=c("re","cr"),k=c(nlevels(pbc2$id),8)),
  sigma ~ 1,
  alpha ~ s(years,k=20),
  dalpha ~ -1
)
## Set the seed for reproducibility.
set.seed(123)

## Estimate model.
# b <- bamlss(f, data = pbc2, family = "jm",
#             timevar = "year", idvar = "id")
b <- readRDS(paste0("/home/alex/Documents/joint_models/JointModel/PBC_analysis",
                    "/pbc_uni_homepage_b.Rds"))
summ <- summary(b$samples[[1]][, grep("\\.accepted", colnames(b$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))

## Set the seed for reproducibility.
set.seed(123)

## Estimate model.
# b <- bamlss(f, data = pbc2, family = "jm",
#             timevar = "year", idvar = "id")
m <- readRDS(paste0("/home/alex/Documents/joint_models/JointModel/PBC_analysis",
                    "/pbc_uni_homepage_m.Rds"))
summ <- summary(m$samples[[1]][, grep("\\.accepted", colnames(m$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))

# Estimate model but now let it converge.
# b <- bamlss(f, data = pbc2, family = "jm",
#             timevar = "year", idvar = "id", maxit = 1500)
b_conv <- readRDS(paste0("/home/alex/Documents/joint_models/JointModel/",
                         "PBC_analysis/pbc_uni_homepage_b_converged.Rds"))
summ <- summary(b_conv$samples[[1]][, grep("\\.accepted", 
                                           colnames(b_conv$samples[[1]]))])
matrix(summ$statistics[, 1], ncol = 1, 
       dimnames = list(names(summ$statistics[, 1]), "Mean"))


# Multivariate Models with Nu Update --------------------------------------

# Estimate the model using estimated FPCs
few_obs <- apply(table(p_long$id, p_long$marker), 1,
                 function (x) any(x < 2))
long_obs <- p_long %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(obstime), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1 * max(p_long$survtime)) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))
mfpca_est <- JMbamlss:::preproc_MFPCA(p_long %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "y ~ 1 + s(obstime) + sex",
                                      save_uniFPCA = TRUE)
# saveRDS(mfpca_est, file = "~/Downloads/mfpca_est.rds")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long <- JMbamlss:::attach_wfpc(mfpca_est, p_long, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker) + sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
set.seed(1714)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 1500, nu = 1, update_nu = TRUE,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

set.seed(1714)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 997, nu = 1, update_nu = TRUE,
                sampler = FALSE, verbose = TRUE,
                par_trace = TRUE)
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:997) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 18, linetype = "dotted") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")
saveRDS(b_est, file = paste0("~/Documents/joint_models/JointModel/",
                             "PBC_analysis/pbc_mul_nuupdate.Rds"))



# Multivariate Analysis on Log-Scale --------------------------------------


# Estimate the model using estimated FPCs
few_obs <- apply(table(p_long$id, p_long$marker), 1,
                 function (x) any(x < 2))
long_obs <- p_long %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(obstime), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1 * max(p_long$survtime)) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))
mfpca_est <- JMbamlss:::preproc_MFPCA(p_long %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "logy ~ 1 + s(obstime) + sex",
                                      save_uniFPCA = TRUE)
# saveRDS(mfpca_est, file = "~/Downloads/mfpca_est.rds")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- JMbamlss:::attach_wfpc(mfpca_est, p_long, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
set.seed(1604)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 1500, verbose = TRUE, 
                par_trace = TRUE, accthreshold = -1)
# saveRDS(b_est, file = paste0("~/Documents/joint_models/JointModel/",
#                              "PBC_analysis/pbc_mul_nuupdate_log.Rds"))
# saved on clapton in JMbamlss/inst/objects/230615_pbc_est.Rds


# Multivariate Analysis on Log Scale and Further Scaling ------------------

p_long$logy10 <- p_long$logy / sqrt(10)

# Estimate the model using estimated FPCs
few_obs <- apply(table(p_long$id, p_long$marker), 1,
                 function (x) any(x < 2))
long_obs <- p_long %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(obstime), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1 * max(p_long$survtime)) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))
mfpca_est <- JMbamlss:::preproc_MFPCA(p_long %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "logy10 ~ 1 + s(obstime) + sex",
                                      save_uniFPCA = TRUE)
# saveRDS(mfpca_est, file = "~/Downloads/mfpca_est.rds")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.975))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- JMbamlss:::attach_wfpc(mfpca_est, p_long, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = logy10) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker) + sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
set.seed(1604)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 1500, update_nu = TRUE,
                verbose = TRUE)
# Numerical Problems in Sampler
# saveRDS(b_est, file = paste0("~/Documents/joint_models/JointModel/",
#                              "PBC_analysis/pbc_mul_nuupdate_log10.Rds"))


# For Overleaf ------------------------------------------------------------

data("pbc2", package = "JMbayes")
pbc2$years <- pbc2$years + sqrt(.Machine$double.eps)
pbc2$year <- pbc2$year + sqrt(.Machine$double.eps)
## Set up the model formula including
## functional random intercepts using ti().
f <- list(
  Surv2(years, status2, obs = log(serBilir)) ~ s(years,k=20),
  gamma ~ s(age,k=20) + drug + sex,
  mu ~ ti(id,bs="re") + 
    ti(year,k=20) + 
    ti(id,year,bs=c("re","cr"),k=c(nlevels(pbc2$id),8)),
  sigma ~ 1,
  alpha ~ s(years,k=20),
  dalpha ~ -1
)
## Set the seed for reproducibility.
set.seed(123)
## Estimate model.
b <- bamlss(f, data = pbc2, family = "jm",
            timevar = "year", idvar = "id")
saveRDS(b, "~/Downloads/pbc2_afterchanges.rds")

data("pbc2", package = "JMbayes")
# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
# Logarithmic scale for all three markers
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years + sqrt(.Machine$double.eps),
         obstime = year + sqrt(.Machine$double.eps), marker = factor(marker),
         logy = log(y)) %>%
  select(id, survtime, event, sex, drug, age, marker, obstime, y, logy) %>%
  arrange(marker, id, obstime) %>%
  na.omit() %>%
  as.data.frame()

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(dplyr::n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis)

f_fri <- list(
  Surv2(survtime, event, obs = logy) ~ s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  mu ~ sex + ti(id, bs = "re") + 
    ti(obstime) + 
    ti(id, obstime, bs = c("re", "cr"), k = c(nlevels(pbc2$id), 8)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1559)
b_bil <- bamlss(f_fri, family = "jm",
                data = p_long %>% filter(marker == "serBilir") %>%
                  droplevels() %>% as.data.frame(),
                timevar = "obstime", idvar = "id", maxit = 1500, verbose = TRUE)
saveRDS(b_bil, "~/Downloads/pbc2_model2.rds")

f_fri_alph <- list(
  Surv2(survtime, event, obs = logy) ~ s(survtime, k = 20),
  gamma ~ s(age, k = 20) + drug + sex,
  mu ~ ti(id, bs="re") + 
    ti(obstime, k=20) + 
    ti(id, obstime, bs = c("re", "cr"), k = c(nlevels(pbc2$id), 8)),
  sigma ~ 1,
  alpha ~ s(survtime, k = 20),
  dalpha ~ -1
)
set.seed(1559)
b_mod1 <- bamlss(f_fri_alph, family = "jm",
                 data = p_long %>% filter(marker == "serBilir") %>%
                   droplevels() %>% as.data.frame(),
                 timevar = "obstime", idvar = "id", maxit = 1500, verbose = TRUE)
saveRDS(b_mod1, "~/Downloads/pbc2_model1.rds")

b_old <- readRDS(paste0("/home/alex/Documents/joint_models/JointModel/PBC_analysis",
                        "/pbc_uni_homepage_b.Rds"))

ri_ids <- grep("\\(id\\)", colnames(b$samples[[1]]))[1:312]
ri_ids_m1 <- grep("\\(id\\)", colnames(b_mod1$samples[[1]]))[1:304]
ri_ids_m2 <- grep("\\(id\\)", colnames(b_bil$samples[[1]]))[1:304]
ri_tau <- unname(summary(b$samples[[1]][,ri_ids[312] + 1])$statistics[1])
ri_tau_o <- unname(summary(b_old$samples[[1]][,ri_ids[312] + 1])$statistics[1])

all_samps <- data.frame(
  id = rep(rep(1:312, each = 1001), 2),
  it = rep(rep(1:1001, 312), 2),
  ri = c(c(b$samples[[1]][, ri_ids]), c(b_old$samples[[1]][, ri_ids])),
  type = factor(rep(c("After Changes", "Before Changes"), each = 1001*312))
)
all_taus <- data.frame(
  x = c(seq(-4, 4, by = 0.1), seq(-4, 4, by = 0.1)),
  y = c(dnorm(seq(-4, 4, by = 0.1), 0, sqrt(ri_tau)),
        dnorm(seq(-4, 4, by = 0.1), 0, sqrt(ri_tau_o))),
  type = factor(c(rep("After Changes", length(seq(-4, 4, by = 0.1))),
                  rep("Before Changes", length(seq(-4, 4, by = 0.1)))))
)

ggplot(all_samps, aes(x = ri)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  ggtitle("Distribution of All Random Intercept Effect Samples",
          paste0("Univariate JM Bilirubin (bamlss), 312 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed") +
  labs(x = "Random Intercept Parameters", y = "Density")

b_bil_old <- readRDS("~/Downloads/pbc2_model2.rds")


# After changing the scaling
b_bil_old <- readRDS("~/Downloads/pbc2_model2.rds")
ri_ids_mold <- grep("\\(id\\)", colnames(b_bil_old$samples[[1]]))[1:304]
ri_ids_mnew <- grep("\\(id\\)", colnames(b_bil$samples[[1]]))[1:304]
ri_tau_o <- unname(summary(b_bil_old$samples[[1]][,ri_ids_mold[304] +
                                                    1])$statistics[1])
ri_tau_n <- unname(summary(b_bil$samples[[1]][,ri_ids_mnew[304] + 
                                                1])$statistics[1])
all_samps <- data.frame(
  id = rep(rep(1:304, each = 1001), 2),
  it = rep(rep(1:1001, 304), 2),
  ri = c(c(b_bil_old$samples[[1]][, ri_ids_mold]),
         c(b_bil$samples[[1]][, ri_ids_mnew])),
  type = factor(rep(c("Submodel 2", "Now"), each = 1001*304),
                levels = c("Submodel 2", "Now"))
)
all_taus <- data.frame(
  x = c(seq(-4, 4, by = 0.1), seq(-4, 4, by = 0.1)),
  y = c(dnorm(seq(-4, 4, by = 0.1), 0, sqrt(ri_tau_o)),
        dnorm(seq(-4, 4, by = 0.1), 0, sqrt(ri_tau_n))),
  type = factor(c(rep("Submodel 2", length(seq(-4, 4, by = 0.1))),
                  rep("Now", length(seq(-4, 4, by = 0.1)))),
                levels = c("Submodel 2", "Now"))
)

ggplot(all_samps, aes(x = ri)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  facet_wrap(~type, scales = "free") +
  theme_bw() +
  ggtitle("Distribution of All Random Intercept Effect Samples",
          paste0("Univariate JM Bilirubin (bamlss), 304 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed") +
  labs(x = "Random Intercept Parameters", y = "Density")

