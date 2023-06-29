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


# Model Fit ---------------------------------------------------------------

f_est <- list(
  Surv2(survtime, event, obs = logy) ~ -1 + 
    s(survtime, k = 10, bs = "ps", xt = list("scale" = FALSE)),
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
# set.seed(1604)
# b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
#                 timevar = "obstime", maxit = 1500, verbose = TRUE, 
#                 par_trace = TRUE, accthreshold = -1)
# saved on clapton in JMbamlss/inst/objects/230615_pbc_est.Rdata

b_est <- readRDS(file.path(server_wd, "JMbamlss/inst/objects",
                           "230615_pbc_est.Rds"))



# Evaluate Model ----------------------------------------------------------


# Acceptance probabilities
acc(b_est)
# xtable::xtable(acc(b_est))

# Trace plot of updating the alpha parameter
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = seq_len(ncol(alpha))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 190, linetype = "dotted") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Alpha Parameters")

# Variance estimates
fpc1 <- grep("\\(id,fpc\\.1\\)", colnames(b_est$samples[[1]]))[1:304]
fpc3 <- grep("\\(id,fpc\\.3\\)", colnames(b_est$samples[[1]]))[1:304]
all_draws <- data.frame(
  id = rep(rep(seq_len(304), each = 1001), times = 2),
  re = c(b_est$samples[[1]][, fpc1], b_est$samples[[1]][, fpc3]),
  type = factor(rep(c("FPC 1", "FPC 3"), each = 1001*304))
)
fpc1_tau2 <- unname(summary(b_est$samples[[1]][, fpc1[304] + 1])$statistics[1])
fpc3_tau2 <- unname(summary(b_est$samples[[1]][, fpc3[304] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-50, 50, by = 0.5), seq(-20, 20, by = 0.1)),
  y = c(dnorm(seq(-50, 50, by = 0.5), 0, sqrt(fpc1_tau2)),
        dnorm(seq(-20, 20, by = 0.1), 0, sqrt(fpc3_tau2))),
  type = factor(c(rep("FPC 1", length(seq(-50, 50, by = 0.5))),
                  rep("FPC 3", length(seq(-20, 20, by = 0.1)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("Distribution of All FPC Score Samples", 
          paste0("MJM PBC (bamlss), 304 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Scores", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")
  


# Analyze Collinearity ----------------------------------------------------

# Argument coll needed
devtools::load_all(".")

set.seed(1604)
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = p_long,
                timevar = "obstime", maxit = 1500, verbose = TRUE,
                par_trace = TRUE, accthreshold = -1, coll = TRUE, 
                sampler = FALSE)


nmarker <- length(mfpca_est$functions)
start <- max(which(sapply(b_est$model.stats$optimizer$coll, is.null))) + 1
stop <- length(b_est$model.stats$optimizer$coll)

# Use data matrix
col_x <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  X <- matrix(x$X, ncol = nmarker)
  svd(crossprod(X))$d
})
condi_x <- apply(col_x, 2, function (x) x[1] / x[nmarker])
ratio_x <- apply(col_x, 2, function (x) x / sum(x))

# Use Hesse matrix
col_i <- sapply(b_est$model.stats$optimizer$coll[start:stop], function (x) {
  H <- matrix(x$I, ncol = nmarker)
  s <- diag(diag(H)^(-1/2))
  svd(s %*% H %*% s)$d
})
condi_i <- apply(col_i, 2, function (x) x[1] / x[nmarker])
ratio_i <- apply(col_i, 2, function (x) x / sum(x))

condi_dat <- data.frame(
  it = rep(start:stop, times = 2),
  vals = c(condi_x, condi_i),
  type = factor(rep(c(0, 1), each = length(start:stop)), labels = c("XTX", "I_S"))
)
ggplot(condi_dat, aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Condition Numbers for PBC Data") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)


ratio_dat <- data.frame(
  it = rep(rep(start:stop, times = nmarker), 2),
  vals = c(ratio_x, ratio_i),
  type = factor(rep(c(0, 1), each = nmarker * length(start:stop)),
                labels = c("XTX", "I_S")),
  ev = factor(seq_len(nmarker), labels = paste("EV", seq_len(nmarker)))
)
ggplot(ratio_dat, aes(x = it, y = vals, col = type,
                      group = interaction(ev, type))) +
  geom_line() +
  theme_bw() +
  ggtitle("Ratio of Eigenvalues for PBC Data") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)
