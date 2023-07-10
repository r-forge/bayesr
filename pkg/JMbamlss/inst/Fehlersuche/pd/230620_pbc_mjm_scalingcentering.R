
# PCB Data Analysis with Centering + Scaling ------------------------------


# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}
server_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.", 
                    "hu-berlin.de,share=volkmana.hub")


# Always
library(survival)
library(JMbayes2)
library(JMbamlss)
library(tidyverse)

# Convenience function
acc <- function(model){
  summ <- summary(model$samples[[1]][, grep("accepted", 
                                            colnames(model$samples[[1]]))])
  matrix(summ$statistics[, 1], ncol = 1, 
         dimnames = list(names(summ$statistics[, 1]), "Mean"))
}



# Prepare the PBC data ----------------------------------------------------


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

# Which longitudinal observations to use for FPCA
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

# Calculating the scaling factors for the uncentered data
mfpca_est <- JMbamlss:::preproc_MFPCA(p_long %>%
                                        filter(id %in% take) %>% 
                                        droplevels(), 
                                      uni_mean = "logy ~ 1 + s(obstime) + sex",
                                      save_uniFPCA = TRUE)
uni_var <- sapply(attr(mfpca_est, "uniFPCA"), function(x) {
   sum(x$evalues) + x$sigma2
})
uni_var_dat <- data.frame(
  marker = as.factor(c("serBilir", "serChol", "SGOT")),
  scale = uni_var
)

# Centering the data OR scaling the data OR both
p_long_cs <- p_long %>%
  left_join(uni_var_dat, by = "marker") %>%
  group_by(marker) %>% 
  mutate(logy_c = logy - mean(logy),
         logy_s = logy / sqrt(scale),
         logy_cs = logy_c / sqrt(scale)) %>%
  ungroup() %>%
  as.data.frame()



# Calculate the different MFPCs -------------------------------------------

mfpca_c <- JMbamlss:::preproc_MFPCA(
  p_long_cs %>% filter(id %in% take) %>% droplevels(), 
  uni_mean = "logy_c ~ 1 + s(obstime) + sex",
  save_uniFPCA = TRUE)
mfpca_s <- JMbamlss:::preproc_MFPCA(
  p_long_cs %>% filter(id %in% take) %>% droplevels(), 
  uni_mean = "logy_s ~ 1 + s(obstime) + sex",
  save_uniFPCA = TRUE)
# MFPC for scaling + centering is the same as only scaling

nfpc_c <- min(which(
  cumsum(mfpca_c$values)/sum(mfpca_c$values) > 0.95))
nfpc_s <- min(which(
  cumsum(mfpca_s$values)/sum(mfpca_s$values) > 0.95))

mfpca_c_list <- lapply(seq_len(nfpc_c), function (i, mfpca = mfpca_c) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
mfpca_s_list <- lapply(seq_len(nfpc_s), function (i, mfpca = mfpca_s) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})

# Prepare objects for model fit
p_c <- JMbamlss:::attach_wfpc(mfpca_c, p_long_cs, n = nfpc_c)
p_s <- JMbamlss:::attach_wfpc(mfpca_s, p_long_cs, n = nfpc_s)



# Model Fit ---------------------------------------------------------------

f_c <- list(
  Surv2(survtime, event, obs = logy_c) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker +",
    paste0(lapply(seq_len(nfpc_c), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_c_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

f_s <- list(
  Surv2(survtime, event, obs = logy_s) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker +",
    paste0(lapply(seq_len(nfpc_s), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_s_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

f_cs <- list(
  Surv2(survtime, event, obs = logy_cs) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker +",
    paste0(lapply(seq_len(nfpc_s), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_s_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Model fit
set.seed(1023)
b_c <- bamlss(f_c, family = JMbamlss:::mjm_bamlss, data = p_c,
                timevar = "obstime", maxit = 1500, verbose = TRUE,
                par_trace = TRUE, accthreshold = -1, coll = TRUE)
saveRDS(b_c, file = file.path(server_wd, "JMbamlss/inst/objects", 
                              "230620_pbc_c.Rds"))
rm(b_c)

set.seed(1023)
b_s <- bamlss(f_s, family = JMbamlss:::mjm_bamlss, data = p_s,
              timevar = "obstime", maxit = 1500, verbose = TRUE,
              par_trace = TRUE, accthreshold = -1, coll = TRUE)
saveRDS(b_s, file = file.path(server_wd, "JMbamlss/inst/objects", 
                              "230620_pbc_s.Rds"))
rm(b_s)

set.seed(1023)
b_cs <- bamlss(f_cs, family = JMbamlss:::mjm_bamlss, data = p_s,
              timevar = "obstime", maxit = 1500, verbose = TRUE,
              par_trace = TRUE, accthreshold = -1, coll = TRUE)
saveRDS(b_cs, file = file.path(server_wd, "JMbamlss/inst/objects", 
                               "230620_pbc_cs.Rds"))


# Why should scaling change stuff? ----------------------------------------

b_cs <- readRDS(file.path(server_wd, "JMbamlss/inst/objects", 
                          "230620_pbc_cs.Rds"))
alpha <- sapply(b_cs$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = seq_len(ncol(alpha))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")

x_s <- p_long_cs %>%
  group_by(marker, id) %>%
  slice(1) %>%
  ungroup()
tapply(x_s$logy_cs, x_s$marker, mean)

p_as <- p_s %>%
  mutate(age_s = scale(age, center = TRUE, scale = TRUE))

f_cs_as <- list(
  Surv2(survtime, event, obs = logy_cs) ~ -1 + 
    s(survtime, k = 20, bs = "ps", xt = list("scale" = FALSE)),
  gamma ~ 1 + age_s + sex,
  as.formula(paste0(
    "mu ~ -1 + marker + s(obstime, by = marker, xt = list('scale' = FALSE))",
    "+ sex:marker +",
    paste0(lapply(seq_len(nfpc_s), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_s_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
set.seed(1023)
b_cs_as <- bamlss(f_cs_as, family = JMbamlss:::mjm_bamlss, data = p_as,
                  timevar = "obstime", maxit = 1500, verbose = TRUE,
                  par_trace = TRUE, accthreshold = -1, coll = TRUE, 
                  sampler = FALSE)

alpha <- sapply(b_cs_as$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = seq_len(ncol(alpha))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")
# also diverges

start <- max(which(sapply(b_cs$model.stats$optimizer$coll, is.null))) + 1
stop <- length(b_cs$model.stats$optimizer$coll)
full_xlong <-  do.call(rbind, 
                       Map(cbind, it = start:stop, eta_mu = 
                           lapply(b_cs$model.stats$optimizer$coll[start:stop], 
                                  function (x) x$X))) %>%
  cbind(p_s %>% group_by(id) %>% slice(1) %>% ungroup() %>% 
          select(survtime)) %>%
  cbind(p_s %>% group_by(marker, id) %>% slice(1) %>% ungroup() %>% 
          select(marker, sex, age, event))

library(manipulate)
manipulate({
  ggplot(full_xlong %>% filter(it == x), 
         aes(x = survtime, y = eta_mu, col = as.factor(event))) + 
    geom_point() +
    facet_grid(~marker) +
    theme_bw() +
    ylim(c(-5, 3))
  }, x = slider(start, stop))

manipulate({
  ggplot(full_xlong %>% filter(it == x), 
         aes(x = survtime, y = eta_mu, col = sex)) + 
    geom_point() +
    facet_grid(~marker) +
    theme_bw() +
    ylim(c(-5, 3))
}, x = slider(start, stop))

manipulate({
  ggplot(full_xlong %>% filter(it == x), 
         aes(x = survtime, y = eta_mu, col = age)) + 
    geom_point() +
    facet_grid(~marker) +
    theme_bw() +
    ylim(c(-5, 3))
}, x = slider(start, stop))


std_steps <- do.call(rbind, 
        Map(cbind, it = start:stop,
              lapply(b_cs$model.stats$optimizer$coll[start:stop], 
                     function (x) {
                       X <- data.frame(matrix(x$X, ncol = 3))
                       data.frame(mean = colMeans(X),
                                  sd = sapply(X, sd))
                     }))) %>%
  cbind(data.frame(marker = as.factor(c("serBilir", "serChol", "SGOT"))))
tapply(std_steps$mean, std_steps$marker, summary)
tapply(std_steps$sd, std_steps$marker, summary)
