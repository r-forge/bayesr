# This file fits several JMbamlss models with different MFPCs. They differ in 
# the number of basis functions used in the estimation of the MFPCs. Finally,
# use the DIC to select the number of basis functions.

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


# Fit the JMbamlss Model --------------------------------------------------

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

# Loop for different number of basis functions
for (nbas in 4:9) {
  
  # Estimate the MFPCA
  mfpca_est <- JMbamlss:::preproc_MFPCA(p_long %>%
                                          filter(id %in% take) %>% 
                                          droplevels(), 
                                        uni_mean = paste0("logy ~ 1 + s(obstime) +
                                                        sex + s(age)"),
                                        fve_uni = 1, nbasis = nbas,
                                        save_uniFPCA = TRUE)
  
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
  t_est <- system.time(
    bamlss_fit <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                         timevar = "obstime", maxit = 1500,
                         n.iter = 12000L, burnin = 2000L, thin = 5)
  )
  attr(bamlss_fit, "mfpca") <- mfpca_est
  attr(bamlss_fit, "estimation_t") <- t_est
  saveRDS(bamlss_fit, file = file.path(results_wd, "inst", "objects", 
                                       "pbc_analysis", 
                                       paste0("pbc_fit_", nbas, "_99.rds")))
  rm(bamlss_fit)
}

dics <- vector(length = 6)
for (nbas in 4:9) {
  bamlss_fit <- readRDS(file.path(results_wd, "inst", "objects", 
                                  "pbc_analysis", 
                                  paste0("pbc_fit_", nbas, "_99.rds")))
  
  dics[nbas-3] <- bamlss_fit$model.stats$sampler$DIC
}
# > dics
# [1]   47.32648  -23.67283 -412.70238 -479.86343 -249.84627 -263.44913