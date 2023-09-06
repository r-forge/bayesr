# Set up R session --------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if(location == "server_linux") "./simulation"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss")
}



# Always
library(survival)
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
library(JMbamlss)



# Prepare Data ------------------------------------------------------------


# Cox model for the composite event death or transplantation
pbc2$event <- as.numeric(pbc2$status != 'alive')

# Longitudinal format for pbc data
# Logarithmic scale for all three markers
p_long <- pbc2 %>%
  pivot_longer(c(serBilir, serChol, SGOT, albumin), names_to = "marker",
               values_to = "y") %>%
  mutate(survtime = years + sqrt(.Machine$double.eps),
         obstime = year + sqrt(.Machine$double.eps), marker = factor(marker),
         logy = log(y)) %>%
  select(id, survtime, event, sex, drug, age, marker, obstime, y, logy) %>%
  arrange(marker, id, obstime) %>%
  na.omit() %>%
  as.data.frame()

# Cap the survival times at the maximum observation times
max_obstime <- max(p_long$obstime)
p_long <- p_long %>%
  mutate(survtime = ifelse(survtime > max_obstime, max_obstime, survtime))

# Remove patients that do not have at least one serChol
which_mis <- p_long %>%
  filter(marker == "serChol") %>%
  group_by(id) %>%
  summarise(dplyr::n()) %>%
  select(id) %>%
  unlist()
p_long <- p_long %>%
  filter(id %in% which_mis) %>%
  droplevels()

# Data for JMbayes2
# Use logy for logarithmic transform
p_long_jmb <- p_long %>%
  pivot_wider(id_cols = c(id, obstime, survtime, event, sex, age),
              values_from = logy, names_from = marker) 

p_long_id <- p_long_jmb %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()



# Fit the JMbamlss Model --------------------------------------------------


# Model using 10 Basis Functions ------------------------------------------

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
                                      uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                        sex + s(age)"),
                                      fve_uni = 1,
                                      save_uniFPCA = TRUE)
# saveRDS(mfpca_est, file = "~/Downloads/mfpca_est.rds")
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
  cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.99))
mfpca_est_list <- lapply(vals, function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- JMbamlss:::attach_wfpc(mfpca_est, p_long, n = length(vals))

# Model formula
f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + s(age) + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker + s(age, by = marker, xt = list('scale' = FALSE)) + ",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1252)
t_99 <- system.time(
  bamlss_fit <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                       timevar = "obstime", maxit = 1500,
                       n.iter = 12000L, burnin = 2000L, thin = 5)
)
saveRDS(bamlss_fit, file = file.path(results_wd, "inst", "objects", 
                                     "pbc_analysis", "pbc_fit_99.rds"))
rm(bamlss_fit)
# Also: fit for truncation order 99.9% and 100%
# # 99.9
# nfpc999 <- min(which(
#   cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.999))
# f_est999 <- list(
#   Surv2(survtime, event, obs = logy) ~ -1 + 
#     s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
#   gamma ~ 1 + s(age) + sex,
#   as.formula(paste0(
#     "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
#     "+ sex:marker + s(age, by = marker, xt = list('scale' = FALSE)) + ",
#     paste0(lapply(seq_len(nfpc999), function(x) {
#       paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
#              "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
#     }), collapse = " + "))),
#   sigma ~ -1 + marker,
#   alpha ~ -1 + marker
# )
# set.seed(1252)
# t_999 <- system.time(
#   bamlss_fit <- bamlss(f_est999, family = JMbamlss:::mjm_bamlss, data = p_long,
#                        timevar = "obstime", maxit = 1500,
#                        n.iter = 12000L, burnin = 2000L, thin = 5)
# )
# saveRDS(bamlss_fit, file = file.path(results_wd, "inst", "objects", 
#                                      "pbc_analysis", "pbc_fit_999.rds"))
# rm(bamlss_fit)
# 100
nfpc1 <- length(vals)
f_est1 <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + s(age) + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker + s(age, by = marker, xt = list('scale' = FALSE)) + ",
    paste0(lapply(seq_len(nfpc1), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
set.seed(1252)
t_1 <- system.time(
  bamlss_fit <- bamlss(f_est1, family = JMbamlss:::mjm_bamlss, data = p_long,
                       timevar = "obstime", maxit = 1500,
                       n.iter = 12000L, burnin = 2000L, thin = 5)
)
saveRDS(bamlss_fit, file = file.path(results_wd, "inst", "objects", 
                                     "pbc_analysis", "pbc_fit_1.rds"))
rm(bamlss_fit)



# JMbayes2 Model ----------------------------------------------------------

# Get the quantile-based knots for comparability
kn <- mgcv::smoothCon(mgcv::s(survtime, k = 10, bs = "ps"), 
                      data = p_long)[[1]]$knots

set.seed(1205)
# Cox Model
t_jmb <- system.time({
  
  CoxFit <- coxph(Surv(survtime, event) ~ ns(age, df = 3) + sex,
                  data = p_long_id)
  
  # Univariate longitudinal models
  fm1ns <- lme(albumin ~ ns(obstime, df = 3) + sex + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit,
               random = ~ ns(obstime, df = 2) | id)
  fm2ns <- lme(serBilir ~ ns(obstime, df = 3) + sex + ns(age, df = 3), 
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  fm3ns <- lme(serChol ~ ns(obstime, df = 3) + sex + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  fm4ns <- lme(SGOT ~ ns(obstime, df = 3) + sex + ns(age, df = 3),
               data = p_long_jmb, na.action = na.omit, 
               random = ~ ns(obstime, df = 2) | id)
  
  
  # Multivariate Joint Model
  # the joint model that links all sub-models
  jointFit <- jm(CoxFit, list(fm1ns, fm2ns, fm3ns, fm4ns), time_var = "obstime",
                 n_iter = 12000L, n_burnin = 2000L, n_thin = 5L,
                 cores = 1, n_chains = 1, 
                 GK_k = 7, Bsplines_degree = 3, diff = 2, knots = list(kn),
                 save_random_effects = TRUE)
  
})
saveRDS(jointFit, file = file.path(results_wd, "inst", "objects", 
                                   "pbc_analysis", "pbc_fit_jmb.rds"))

# Sensitivity Analysis with weighted MFPCs --------------------------------


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
mfpca_w <- JMbamlss:::preproc_MFPCA(p_long %>% filter(id %in% take) %>% 
                                      droplevels(),
                                    uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                      sex + s(age)"),
                                    fve_uni = 1, nbasis = 7, weights = TRUE,
                                    save_uniFPCA = TRUE)
vals <- which(mfpca_w$values > 0)
nfpc <- min(which(
  cumsum(mfpca_w$values[vals])/sum(mfpca_w$values[vals]) > 0.99))
mfpca_w_list <- lapply(vals, function (i, mfpca = mfpca_w) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_long[, grepl("fpc", colnames(p_long))] <- NULL
p_long <- JMbamlss:::attach_wfpc(mfpca_w, p_long, n = length(vals))

# Model formula
f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + s(age) + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker + s(age, by = marker, xt = list('scale' = FALSE)) + ",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_w_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1252)
t_w <- system.time(
  bamlss_wei <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                       timevar = "obstime", maxit = 1500,
                       n.iter = 12000L, burnin = 2000L, thin = 5)
)
saveRDS(bamlss_wei, file = file.path(results_wd, "inst", "objects", 
                                     "pbc_analysis", "pbc_fit_w.rds"))


# Computation Times -------------------------------------------------------

# bamlss_1
# > t_1
# user    system   elapsed 
# 11983.029    15.524 11995.517 
# After changes (while working)
# > t_1
# user    system   elapsed 
# 14762.260    54.023 14812.989 


# bamlss_999
# > t_999
# user   system  elapsed 
# 10376.79     7.12 10381.26 

# bamlss_99
# > t_99
# user   system  elapsed 
# 7142.438    3.056 7143.590 
# After changes (while working)
# > t_99
# user   system  elapsed 
# 7534.451   35.993 7568.498 

# jmb
# > t_jmb
# user  system elapsed 
# 723.656   2.784 726.288 
# After changes
# > t_jmb
# user  system elapsed 
# 197.138   0.848 197.960 

# Weighted (while working)
# > t_w
# user   system  elapsed 
# 9297.963  141.809 9436.353 

# Plot Observations -------------------------------------------------------

# Create data set for prediction
id_grid <- sapply(p_long_id$survtime, function (x) {
  c(seq(0 + sqrt(.Machine$double.eps), x, by = 0.2), x)
})
id_grid_length <- sapply(id_grid, length)
newdat <- p_long_id[rep(seq_len(nrow(p_long_id)), id_grid_length), ] %>%
  mutate(obstime = unlist(id_grid)) %>%
  pivot_longer(cols = c(albumin, serBilir, serChol, SGOT), 
               names_to = "marker") %>%
  mutate(marker = as.factor(marker),
         value = NULL) %>%
  arrange(marker, id, obstime) %>%
  droplevels() %>%
  as.data.frame()

# Prediction of bamlss model fit
preds <- cbind(
  newdat,
  predict(bamlss_fit, newdata = newdat, model = "mu", FUN = c95))

# Prediction of weighted model fit
preds_w <- predict(bamlss_wei, newdata = newdat, model = "mu", FUN = c95)
colnames(preds_w) <- paste0("w", colnames(preds_w))
preds <- cbind(preds, preds_w)

# Prediction of JMbayes2 model fit
# ns() objects have to be predicted from original data set otherwise the values
# do not align
X <- lapply(split(newdat, newdat$marker), function(x) {
  dat <- cbind(predict(ns(p_long_jmb$obstime, df = 3), x$obstime), 
               model.matrix(~ sex, x),
               predict(ns(p_long_jmb$age, df = 3), x$age))
  colnames(dat)[1:3] <- paste0("ns(obstime, df = 3)", 1:3)
  colnames(dat)[6:8] <- paste0("ns(age, df = 3)", 1:3)
  dat[, c(4, 1:3, 5:8)]
})
Z <- lapply(split(newdat, newdat$marker), function(x) {
  dat <- cbind(model.matrix(~ 1, x),
               predict(ns(p_long_jmb$obstime, df = 2), x$obstime))
  colnames(dat)[2:3] <- paste0("ns(obstime, df = 2)", 1:2)
  dat
})
B <- jointFit$mcmc$b[[1]]
n_re <- ncol(Z[[1]])
mcmc_mu <- do.call(rbind, lapply(seq_len(4), function (dim) {
  idL <- unclass(droplevels(split(newdat$id, newdat$marker)[[dim]]))
  tcrossprod(X[[dim]], jointFit$mcmc[[paste0("betas", dim)]][[1]]) +
    t(sapply(seq_len(nrow(Z[[dim]])), function (i) {
      Z[[dim]][i, ] %*% B[idL[i], (dim - 1)*n_re + seq_len(n_re), ]
    }))
}))
jmb_preds <- data.frame("jmb_2.5" = apply(mcmc_mu, 1, quantile, probs = 0.025),
                        "jmb_Mean" = rowMeans(mcmc_mu),
                        "jmb_97.5" = apply(mcmc_mu, 1, quantile, probs = 0.975))
preds <- cbind(preds, jmb_preds)

# Reform the data set for nicer plotting
preds_plot <- pivot_longer(preds, cols = c("Mean", "wMean", "jmb_Mean"), 
                           names_to = "Model", values_to = "fit") %>%
  mutate(Model = factor(Model, levels = c("Mean", "wMean", "jmb_Mean"), 
                        labels = c("bamlss", "bamlss (w)", "JMB")))

# Best of two:
# set.seed(833)
# set.seed(1106)
# ids <- sample(levels(p_long$id), size = 5)
# ids <- c("40", "54", "114", "183")
ids <- c(40, 54, 114, 204)

ggplot(preds_plot %>% 
         filter(id %in% ids) %>%
         droplevels() %>%
         mutate(id = factor(id, labels = paste("Subject", levels(id)))),
       aes(x = obstime, y = fit)) + 
  geom_line(aes(linetype = Model, color = Model)) +
  theme_bw() +
  facet_grid(marker ~ id, scales = "free_y") +
  geom_point(data = p_long %>% 
               filter(id %in% ids) %>%
               droplevels() %>%
               mutate(id = factor(id, labels = paste("Subject", levels(id)))),
             aes(y = logy), alpha = 0.5) +
  labs(y = NULL, x = "Time")
# save as 4x8

# library(manipulate)
# manipulate(ggplot(preds_plot %>% filter(id %in% ((ids-1)*4 + 1:4)),
#                   aes(x = obstime, y = fit)) +
#              geom_line(aes(linetype = Model, color = Model)) +
#              theme_bw() +
#              facet_grid(marker ~ id, scales = "free_y") +
#              geom_point(data = p_long %>% filter(id %in% ((ids-1)*4 + 1:4)),
#                         aes(y = logy), alpha = 0.5) +
#              labs(y = NULL, x = "Time"),
#            ids = slider(1, 76))


# Plot Estimated Eigenfunctions -------------------------------------------

# Estimate the MFPCA
mfpca_4 <- JMbamlss:::preproc_MFPCA(p_long %>% filter(id %in% take) %>% 
                                      droplevels(),
                                    uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                      sex + s(age)"),
                                    fve_uni = 1, nbasis = 4,
                                    save_uniFPCA = TRUE)
mfpca_7 <- JMbamlss:::preproc_MFPCA(p_long %>% filter(id %in% take) %>% 
                                      droplevels(),
                                    uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                      sex + s(age)"),
                                    fve_uni = 1, nbasis = 7,
                                    save_uniFPCA = TRUE)
mfpca_1 <- JMbamlss:::preproc_MFPCA(p_long %>% filter(id %in% take) %>% 
                                      droplevels(),
                                    uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                      sex + s(age)"),
                                    fve_uni = 1, nbasis = 10,
                                    save_uniFPCA = TRUE)
mfpca_w <- JMbamlss:::preproc_MFPCA(p_long %>% filter(id %in% take) %>% 
                                      droplevels(),
                                    uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                      sex + s(age)"),
                                    fve_uni = 1, nbasis = 7, weights = TRUE,
                                    save_uniFPCA = TRUE)


# Extract the information of first two Eigenfunctions and convert to data frame
refObj <- extractObs(mfpca_7$functions, 1:7)
fpcs <- lapply(list(mfpca_4, mfpca_7, mfpca_1, mfpca_w), function(x) {
  dat <- extractObs(x$functions, 1:7)
  dat <- flipFuns(refObj, dat)
  multifammPaper::funData2DataFrame(dat)
})
fpc_dat <- do.call(
  rbind, Map(cbind,
             nbasis = factor(1:4, labels = c("4B", "7B", "10B", "7Bw")),
             fpcs)) %>%
  mutate(dim = factor(dim,
                      labels = c("albumin", "serBilir", "serChol", "SGOT")),
         obs = factor(obs, labels = paste("MFPC", 1:7)))

ggplot(fpc_dat, aes(x = t, y = y, color = nbasis)) +
  geom_line(aes(linetype = nbasis)) +
  facet_grid(obs ~ dim, scales = "free_y") +
  theme_bw() +
  labs(y = NULL, x = "Time", color = NULL, linetype = NULL)
# 6x12


# Description of MFPCA ----------------------------------------------------

uniFPCA <- attr(mfpca_1, "uniFPCA")
sapply(uniFPCA, "[[", "npc")



# Table for Estimated Parameters ------------------------------------------

jointFit$statistics$Mean$gammas["sexfemale"]
aj_mean <- jointFit$statistics$Mean$alphas
aj_low <- jointFit$statistics$CI_low$alphas
aj_up <- jointFit$statistics$CI_upp$alphas

alpha_b <- summary(samples(bamlss_fit)[, grep("alpha\\.",
                                              colnames(samples(bamlss_fit)))])
alpha_bw <- summary(samples(bamlss_wei)[, grep("alpha\\.",
                                               colnames(samples(bamlss_wei)))])
ab_mean <- alpha_b$statistics[1:4, 1]
ab_low <- alpha_b$quantiles[1:4, 1]
ab_up <- alpha_b$quantiles[1:4, 3]
abw_mean <- alpha_bw$statistics[1:4, 1]
abw_low <- alpha_bw$quantiles[1:4, 1]
abw_up <- alpha_bw$quantiles[1:4, 3]

cat(mapply(function(marker, a_mean, a_low, a_up, b_mean, b_low, b_up, c_mean,
                    c_low, c_up){
  paste(marker, 
        "&$", round(a_mean, digits = 2), "$ $[", round(a_low, digits = 2), 
            "$; $", round(a_up, digits = 2), "]$ & ",
        "$", round(b_mean, digits = 2), "$ $[", round(b_low, digits = 2),
            "$; $", round(b_up, digits = 2), "]$ & ",
        "$", round(c_mean, digits = 2), "$ $[", round(c_low, digits = 2),
            "$; $", round(c_up, digits = 2), "]$ \\\\ ")
}, c("albumin", "serBilir", "serChol", "SGOT"), ab_mean, ab_low, ab_up, 
abw_mean, abw_low, abw_up, aj_mean, aj_low, aj_up))
