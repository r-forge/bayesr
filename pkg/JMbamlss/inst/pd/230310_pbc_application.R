
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
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years, obstime = year, marker = factor(marker)) %>%
  select(id, survtime, event, sex, age, marker, obstime, y) %>%
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
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

  
# Standardize the data
univars <- data.frame(marker = factor(c("serBilir", "serChol", "SGOT")),
                      ev_std = sqrt(sapply(attr(mfpca_est, "uniFPCA"), 
                                           function(x) sum(x$evalues))))
p_long_std <- p_long %>%
  left_join(univars, by = "marker") %>%
  mutate(y = y/ev_std)
mfpca_std <- JMbamlss:::preproc_MFPCA(p_long_std %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "y ~ 1 + s(obstime) + sex",
                                      save_uniFPCA = TRUE)
p_long_std[, grep("fpc\\.", names(p_long_std))] <- NULL
vals <- which(mfpca_std$values > 0)
nfpc <- min(which(
  cumsum(mfpca_std$values[vals])/sum(mfpca_std$values[vals]) > 0.975))
p_long_std <- JMbamlss:::attach_wfpc(mfpca_std, p_long_std, n = nfpc)

mfpca_std_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_std) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_std <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker) + sex:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_std_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

b_std <- bamlss(f_std, family = JMbamlss:::mjm_bamlss, data = p_long_std,
                timevar = "obstime", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)



# Univariate Models -------------------------------------------------------

# Random intercept + random slope
f_uni <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + age + sex,
  mu ~ sex + s(obstime) + s(id, bs = "re") + s(obstime, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_bil <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "serBilir"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_cho <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "serChol"),
                timevar = "obstime", idvar = "id", maxit = 5000,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)
b_sgo <- bamlss(f_uni, family = "jm",
                data = p_long %>% filter(marker == "SGOT"),
                timevar = "obstime", idvar = "id", maxit = 1500,
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE)

summary(b_bil$samples[[1]][,grep("\\.alpha", colnames(b_cho$samples[[1]]))])
summary(b_cho$samples[[1]][,grep("\\.alpha", colnames(b_cho$samples[[1]]))])
summary(b_sgo$samples[[1]][,grep("\\.alpha", colnames(b_cho$samples[[1]]))])
# all three unproblematic
