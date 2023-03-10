# Set up R session --------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}
server_wd <- switch(location,
                    "workstation" = paste0("/run/user/1000/gvfs/smb-share:",
                                           "server=clapton.wiwi.hu-berlin.de,",
                                           "share=volkmana.hub/JMbamlss/",
                                           "simulation/"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation/",
                    "server_windows" = "H:/JMbamlss/simulation/")


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

i <- 1
set.seed(i)
setting <- "scen_mauff/"


# Load the data
load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))

# Analysis of the data
ggplot(dat.id, aes(x = Time, colour = as.factor(event))) + 
  geom_histogram(aes(y = after_stat(density), fill = as.factor(event)), 
                 alpha = 0.3, position = "identity") + 
  geom_density() +
  theme_bw() +
  labs(y = "Density", colour = "Event", fill = "Event") +
  ggtitle("Distribution of Event Times in Mauff Simulation",
          "Simulated Data Set 1")
ggplot(dat, aes(x = year, colour = as.factor(event))) + 
  geom_histogram(aes(y = after_stat(density), fill = as.factor(event)), 
                 alpha = 0.3, position = "identity") + 
  geom_density() +
  theme_bw() +
  labs(y = "Density", colour = "Event", fill = "Event") +
  ggtitle("Distribution of Observation Times in Mauff Simulation",
          "Simulated Data Set 1")
helpdat <- dat %>% 
  pivot_longer(y1:y6, names_to = "marker", values_to = "y") %>%
  mutate(marker = factor(marker, labels = paste0("m", 1:6)),
         id = factor(id))
ggplot(helpdat,
       aes(x = year, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories",
          "Simulated Data Set 1")




# Univariate JM -----------------------------------------------------------

dat_uni <- helpdat %>%
  filter(marker == "m1")

f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re") + s(year, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_uni <- bamlss(f_uni, family = "jm", data = dat_uni,
                timevar = "year", idvar = "id", maxit = 1500, n.iter = 900,
                burnin = 1000, thin = 3, verbose = TRUE)


f_uni_fre <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re") + s(id, year, bs = c("re", "cr"),
                                    k = c(500, 8)),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_uni_fre <- bamlss(f_uni_fre, family = "jm", data = dat_uni,
                    timevar = "year", idvar = "id", maxit = 1500, n.iter = 900,
                    burnin = 1000, thin = 3, verbose = TRUE)

# Univariate model with new data set
d_rirs <- readRDS(paste0(server_wd, setting, "/data/data", i, ".rds"))
f_uni <- list(
  Surv2(Time1cens, event, obs = y) ~ -1 + s(Time1cens, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year1cens + s(id, bs = "re") + s(year1cens, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
b_uni <- bamlss(f_uni, family = "jm", data = d_rirs %>% filter(marker == "m1"),
                timevar = "year1cens", idvar = "id", maxit = 1500, n.iter = 900,
                burnin = 1000, thin = 3, verbose = TRUE)


f_mjm <- list(
  Surv2(Time1cens, event, obs = y) ~ -1 + s(Time1cens, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~ year1cens + s(id, bs = "re") + s(id, year1cens, bs = "re"),
  sigma ~ 1,
  alpha ~ 1
)
b_mjm <- bamlss(f_mjm, family = JMbamlss:::mjm_bamlss,
                data = d_rirs %>% filter(marker == "m1"),
                timevar = "year1cens", maxit = 1, sampler = FALSE,
                burnin = 1000, thin = 3, verbose = TRUE, uni = TRUE)
b_MJM <- bamlss(f_mjm, family = JMbamlss:::mjm_bamlss, start = ps,
                data = d_rirs %>% filter(marker == "m1"),
                timevar = "year1cens", maxit = 600, sampler = FALSE,
                burnin = 1000, thin = 3, verbose = TRUE, uni = TRUE)


# Use True MFPC Basis -----------------------------------------------------
f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re") + s(year, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)




# JMbamlss with Weighted MFPCA --------------------------------------------

d_rirs <- readRDS(paste0(server_wd, setting, "/data/data", i, ".rds"))


# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$id, d_rirs$marker), 1, 
                 function (x) any(x < 2))
long_obs <- d_rirs %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(year), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))

mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      time = "year1cens", weights = TRUE,
                                      uni_mean = "y ~ 1 + year1cens",
                                      npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)

uni_norms <- lapply(mfpca_est$functions, norm)
nfpc <- length(vals)
# nfpc <- max(sapply(uni_norms, function (n) {
#   min(which(
#     cumsum(mfpca_est$values[vals]*n)/sum(mfpca_est$values[vals]*n) > 0.9
#   ))
# }))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], 
                         function (i, mfpca = mfpca_est) {
                           list(functions = extractObs(mfpca$functions, i),
                                values = mfpca$values[i])
                         })

# Prepare objects for model fit
d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs, n = nfpc,
                                     obstime = "year1cens")
f_est <- list(
  Surv2(Time1cens, event, obs = y) ~ -1 + s(Time1cens, k = 10, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year1cens:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
#t_est <- system.time(
set.seed(1415)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                timevar = "year1cens", maxit = 0, sampler = FALSE,
                burnin = 1000, thin = 3, verbose = FALSE)
ps <- parameters(b_est)
taus <- grep("tau21", names(ps))
ps[taus] <- ps[taus] /27^2
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                  timevar = "year1cens", maxit = 1500, n.iter = 900,
                  burnin = 1000, thin = 3, start = ps,  verbose = TRUE)
#)
attr(b_est, "comp_time") <- t_est
attr(b_est, "FPCs") <- mfpca_est
attr(b_est, "nfpc") <- nfpc


# -------------- FOLLOWING DOES NOT WORK ----------------------------------
# JMbamlss Model ----------------------------------------------------------

d_rirs <- pivot_longer(dat, y1:y6, names_to = "marker", values_to = "y") %>%
  mutate(marker = factor(marker, labels = paste0("m", 1:6)),
         id = factor(id),
         Time = ifelse(Time/(max(year)) > 1, max(year), Time/(max(year))),
         year = year/max(year)) %>%
  arrange(marker, id, year) %>%
  as.data.frame()

  
# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$id, d_rirs$marker), 1, 
                 function (x) any(x < 2))
long_obs <- d_rirs %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(year), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))

mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      time = "year", method = "PACE",
                                      uni_mean = "y ~ 1 + year + group + year:group",
                                      npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)

uni_norms <- lapply(mfpca_est$functions, norm)
nfpc <- max(sapply(uni_norms, function (n) {
  min(which(
    cumsum(mfpca_est$values[vals]*n)/sum(mfpca_est$values[vals]*n) > 0.9
  ))
}))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], 
                         function (i, mfpca = mfpca_est) {
                           list(functions = extractObs(mfpca$functions, i),
                                values = mfpca$values[i])
                         })

# Prepare objects for model fit
d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs, n = nfpc,
                                     obstime = "year")
f_est <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est <- system.time(
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est, 
                  timevar = "year", maxit = 1500, n.iter = 900,
                  burnin = 1000, thin = 3, verbose = TRUE)
)
attr(b_est, "comp_time") <- t_est
attr(b_est, "FPCs") <- mfpca_est
attr(b_est, "nfpc") <- nfpc



# JMbamlss Model on Shorter Interval --------------------------------------

d_rirs_short <- pivot_longer(dat, y1:y6, names_to = "marker", 
                             values_to = "y") %>%
  mutate(marker = factor(marker, labels = paste0("m", 1:6)),
         id = factor(id),
         Time = ifelse(Time/20 > 1, 1, Time/20),
         year = year/20) %>%
  arrange(marker, id, year) %>%
  as.data.frame()


# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs_short$id, d_rirs_short$marker), 1, 
                 function (x) any(x < 2))
long_obs <- d_rirs_short %>%
  group_by(id, marker) %>%
  summarize(maxobs = max(year), .groups = "drop_last") %>%
  ungroup(marker) %>%
  summarize(minmaxobs = min(maxobs), .groups = "drop_last") %>%
  ungroup() %>%
  filter(minmaxobs > 0.1) %>%
  select(id) %>%
  unlist() %>%
  paste()
take <- intersect(long_obs, paste(which(!few_obs)))

mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs_short %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      time = "year", method = "PACE",
                                      uni_mean = "y ~ 1 + year + group + year:group",
                                      npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)

uni_norms <- lapply(mfpca_est$functions, norm)
nfpc <- max(sapply(uni_norms, function (n) {
  min(which(
    cumsum(mfpca_est$values[vals]*n)/sum(mfpca_est$values[vals]*n) > 0.9
  ))
}))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], 
                         function (i, mfpca = mfpca_est) {
                           list(functions = extractObs(mfpca$functions, i),
                                values = mfpca$values[i])
                         })

# Prepare objects for model fit
d_rirs_short_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs_short, n = nfpc,
                                     obstime = "year")
f_est <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
t_est <- system.time(
  b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_short_est, 
                  timevar = "year", maxit = 1500, n.iter = 900,
                  burnin = 1000, thin = 3, verbose = TRUE)
)
