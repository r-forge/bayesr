

# Collinearity in Mauff ---------------------------------------------------


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
dat_setting <- "scen_mauff"


# Always
library(survival)
library(JMbayes2)
library(MFPCA)
library(tidyverse)
devtools::load_all(".")



# Calculate Iteration 1 with True FPCs ------------------------------------

# Load data
d_rirs <- readRDS(file.path(server_wd, "JMbamlss/simulation", dat_setting,
                            "data/data1.rds"))

# True FPCs
D_num <- matrix(NA, nrow = 12, ncol = 12)
D_num[1:4, 1:4] <- matrix(c(0.919,  0.003,  0.220, -0.008,
                            0.003,  0.006,  0.005, -0.002,
                            0.220,  0.005,  0.421, -0.013,
                            -0.008, -0.002, -0.013,  0.008), byrow = TRUE)
D_num[5:8, 5:8] <- matrix(c(0.551,  0.007, -0.141,  0.015,
                            0.007,  0.014, -0.005, -0.005,
                            -0.141, -0.005,  0.204, -0.042,
                            0.015, -0.005, -0.042,  0.043), byrow = TRUE)  
D_num[9:12, 9:12] <- matrix(c(0.110, -0.028,  0.176,  0.021,
                              -0.028,  0.035, -0.026, -0.003,
                              0.176, -0.026,  3.580,  0.040,
                              0.021, -0.003,  0.040,  0.197), byrow = TRUE)  
D_num[1:4, 5:8] <- matrix(c(0.308,  0.005,  0.265, -0.013,
                            0.012,  0.001,  0.012,  0.002,
                            -0.073, -0.007, -0.172,  0.018,
                            0.007,  0.005,  0.012, -0.013), byrow = TRUE)
D_num[1:4, 9:12] <- matrix(c(0.049,  0.005,  0.124, -0.013,
                             -0.007, -0.005, -0.010,  0.012,
                             0.698,  0.006,  0.680, -0.027,
                             0.056,  0.006,  0.034,  0.004), byrow = TRUE)  
D_num[5:8, 9:12] <- matrix(c(0.095,  0.004, -0.144,  0.032,
                             -0.013,  0.005,  0.035, -0.037,
                             0.826,  0.018, -0.263,  0.024,
                             0.077,  0.019, -0.015,  0.006), byrow = TRUE)
D_num[5:8, 1:4] <- t(D_num[1:4, 5:8])
D_num[9:12, 1:4] <- t(D_num[1:4, 9:12])
D_num[9:12, 5:8] <- t(D_num[5:8, 9:12])

D_numsym <- Matrix::forceSymmetric(D_num) 
D <- D_numsym

# Basis functions on each dimension
seq1 <- seq(0, max(d_rirs$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs)
mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })

# Prepare objects for model fit
nfpc <- 6
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs, n = nfpc,
                                     obstime = "year")

# Reduce nfpc to 6 for faster computation
f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps", 
                                       xt = list("scale = FALSE")),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# # Cap the survival times at the maximum observation times
# max_obstime <- max(d_rirs_tru$year)
# d_rirs_tru <- d_rirs_tru %>%
#   mutate(Time = ifelse(Time > max_obstime, max_obstime, Time))


# Use it with new add-on that Score and Hesse are also returned
set.seed(1432)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "year", verbose = TRUE, sampler = FALSE,
                maxit = 800, std_surv = TRUE, par_trace = TRUE)
set.seed(1213)
b_new <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru,
                timevar = "year", verbose = TRUE, optimizer = FALSE,
                start = parameters(b_est))

# Full model (nfpc = 12) crashes at iteration 541 (alpha starts at 219)


# Analyze Collinearity Model Truncated ------------------------------------

nmarker <- length(mfpca_tru$functions)
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
  ggtitle("Condition Numbers for Mauff It 1 Truncated") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
ggsave("mjm_mauff_trunc_coll_condi.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)


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
  ggtitle("Ratio of Eigenvalues for Mauff It 1 Truncated") +
  labs(x = "Iteration", y = "Eigenval / sum(Eigenval)", col = NULL,
       ev = NULL)
ggsave("mjm_mauff_trunc_coll_ratio.pdf", device = "pdf", path = "plots",
       width = 8, height = 4)


# Further Cutoff (5 FPCs) -------------------------------------------------

nfpc <- 5
f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps", 
                                       xt = list("scale = FALSE")),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "year", maxit = 1500, sampler = FALSE,
                verbose = TRUE, coll = TRUE, par_trace = TRUE)
start <- max(which(sapply(b_est$model.stats$optimizer$coll, is.null))) + 1
stop <- length(b_est$model.stats$optimizer$coll)


# Trace plot of updating the alpha parameter
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
p1 <- ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = seq_len(ncol(alpha))) %>%
         pivot_longer(cols = -it) %>% filter(it %in% 200:225),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(0.1, -0.6, -0.1, -1.41, -1.81, 0.75), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Zoom2: Alpha Parameters Mauff Simulation (Truncated at 5 FPCs)")
p2 <- ggplot(data = data.frame(t(alpha)) %>%
               mutate(it = seq_len(ncol(alpha))) %>%
               pivot_longer(cols = -it) %>% filter(it %in% 200:300),
             aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(0.1, -0.6, -0.1, -1.41, -1.81, 0.75), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Zoom1: Alpha Parameters Mauff Simulation (Truncated at 5 FPCs)")
p3 <- ggplot(data = data.frame(t(alpha)) %>%
               mutate(it = seq_len(ncol(alpha))) %>%
               pivot_longer(cols = -it),
             aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(0.1, -0.6, -0.1, -1.41, -1.81, 0.75), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Alpha Parameters Mauff Simulation (Truncated at 5 FPCs)")
p4 <- gridExtra::grid.arrange(p3, p2, p1, ncol = 1)
ggsave("mjm_mauff_trunc5_alphapar.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 9, plot = p4)


# Trace plot of updating the gamma parameter
gamma <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$gamma$p)
p1 <- ggplot(data = data.frame(t(gamma)) %>%
               mutate(it = seq_len(ncol(gamma))) %>%
               pivot_longer(cols = -it) %>% filter(it %in% 200:225),
             aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(-5.8, 0.5), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Zoom2: Gamma Parameters Mauff Simulation (Truncated at 5 FPCs)")
p2 <- ggplot(data = data.frame(t(gamma)) %>%
         mutate(it = seq_len(ncol(gamma))) %>%
         pivot_longer(cols = -it) %>% filter(it %in% 200:300),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(-5.8, 0.5), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Zoom1: Gamma Parameters Mauff Simulation (Truncated at 5 FPCs)")
p3 <- ggplot(data = data.frame(t(gamma)) %>%
               mutate(it = seq_len(ncol(gamma))) %>%
               pivot_longer(cols = -it),
             aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(-5.8, 0.5), 
             linetype = "dotted") +
  geom_vline(xintercept = 204, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Gamma Parameters Mauff Simulation (Truncated at 5 FPCs)")
p4 <- gridExtra::grid.arrange(p3, p2, p1, ncol = 1)
ggsave("mjm_mauff_trunc5_gammapar.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 9, plot = p4)

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
  type = factor(rep(c(0, 1), each = length(start:stop)), 
                labels = c("XTX", "I_S"))
)
p1 <- ggplot(condi_dat, aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Condition Numbers for Mauff It 1 (Truncated at 5 FPCs)") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
p2 <- ggplot(condi_dat %>% filter(it %in% c(start:300)),
             aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Zoom1: Condition Numbers for Mauff It 1 (Truncated at 5 FPCs)") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
p3 <-ggplot(condi_dat %>% filter(it %in% c(start:225)),
            aes(x = it, y = vals, col = type)) +
  geom_line() +
  theme_bw() +
  ggtitle("Condition Numbers for Mauff It 1 (Truncated at 5 FPCs)") +
  labs(x = "Iteration", y = "max(Eigenval) / min(Eigenval)", col = NULL)
p4 <- gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
ggsave("mjm_mauff_trunc5_coll_condi.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 9, plot = p4)


# Trace plot of updating the lambda part
lambda <- sapply(b_est$model.stats$optimizer$par_trace, 
                 function(x) x$lambda$s[[1]])
ggplot(data = data.frame(t(lambda)) %>%
         select(-tau21, -edf) %>%
         mutate(it = seq_len(ncol(lambda))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = start, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  theme(legend.position = "None") +
  ggtitle("Lambda Parameters Mauff Simulation (Truncated at 5 FPCs)")
ggsave("mjm_mauff_trunc5_lambdapar.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 3)

# Trace plot of updating the mu parameter part
mu <- sapply(b_est$model.stats$optimizer$par_trace, 
             function(x) x$mu$p)
ggplot(data = data.frame(t(mu)) %>%
         mutate(it = seq_len(ncol(gamma))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = start, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Mu Parameters Mauff Simulation (Truncated at 5 FPCs)")
ggsave("mjm_mauff_trunc5_muppar.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 3)

# Longitudinal fixed effects
mus1 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[1]])
ggplot(data = data.frame(t(mus1)) %>%
         select(-tau21, -edf) %>%
         mutate(it = seq_len(ncol(lambda))) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = start, linetype = "dashed") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  theme(legend.position = "None") +
  ggtitle("Mu FPC1 Mauff Simulation (Truncated at 5 FPCs)")
ggsave("mjm_mauff_trunc5_mufpc1par.pdf", device = "pdf", path = "~/Downloads",
       width = 8, height = 3)
