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
d_rirs <- readRDS(paste0(server_wd, setting, "/data/data", i, ".rds"))



# True FPCs ---------------------------------------------------------------

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

D_numsym <- forceSymmetric(D_num) 
D <- D_numsym

# Basis functions on each dimension
seq1 <- seq(0, max(d_rirs$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs)



# TRUE MFPC Model ---------------------------------------------------------

mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })

# Prepare objects for model fit
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs, n = 12,
                                     obstime = "year")
# Reduce nfpc to 5 for faster computation
f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps", 
                                       xt = list("scale = FALSE")),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_along(1:5), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]], 'scale' = FALSE))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

# Use it with new add-on that Score and Hesse are also returned
set.seed(i)
b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "year", maxit = 1500,
                #burnin = 1000, thin = 3,
                verbose = TRUE, sampler = FALSE,
                par_trace = TRUE)
# Doesn't work, also for fewer MFPCs
saveRDS(b_est, file = paste0(server_wd, setting, "/bamlss_est/b_trunc", i,
                             ".rds"))


# Parameter Traceplots ----------------------------------------------------


# Parameter values
lambda <- sapply(b_est$model.stats$optimizer$par_trace, 
                 function(x) x$lambda$s[[1]])
gamma <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$gamma$p)
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
muf1 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[1]])
muf2 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[2]])
muf3 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[3]])
muf4 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[4]])
muf5 <- sapply(b_est$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[5]])
mup <- sapply(b_est$model.stats$optimizer$par_trace, 
              function(x) x$mu$p)
sigma <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$sigma$p)

ggplot(data = data.frame(t(lambda)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Parameters")

ggplot(data = data.frame(t(gamma)) %>%
         mutate(it = 1:79) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Gamma Parameters")

ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:1500) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")

ggplot(data = data.frame(t(muf1)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu FPC1 Parameters")

ggplot(data = data.frame(t(muf2)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu FPC2 Parameters")

ggplot(data = data.frame(t(muf3)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu FPC3 Parameters")

ggplot(data = data.frame(t(muf4)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu FPC4 Parameters")

ggplot(data = data.frame(t(muf5)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu FPC5 Parameters")

ggplot(data = data.frame(t(mup)) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Parameters")

ggplot(data = data.frame(t(sigma)) %>%
         mutate(it = 1:80) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Sigma Parameters")



# Score Trace Plots -------------------------------------------------------
library(gridExtra)

lambda_l <- lapply(b_est$model.stats$optimizer$updates$lambda[[2]], 
                   function(x) {
 data.frame(
   "par" = paste0("par", seq_along(x$xscore)),
   "s_like" = x$xscore0, 
   "s_prio" = x$xscore - x$xscore0, 
   "s_post" = x$xscore) 
})
lambda <- do.call(rbind, Map(cbind, it = seq_along(lambda_l), lambda_l))

ggplot(data = lambda,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Score")

ggplot(data = lambda,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-20, 20)) +
  ggtitle("Lambda Score")

ggplot(data = lambda,
       aes(x = it, y = s_prio, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Score")

p_score_l <- ggplot(data = lambda,
       aes(x = it, y = s_post, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-20, 20)) +
  ggtitle("Lambda Score")

##---
gamma_l <- lapply(b_est$model.stats$optimizer$updates$gamma[[1]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
gamma <- do.call(rbind, Map(cbind, it = seq_along(gamma_l), gamma_l))
p_score_g <- ggplot(data = gamma,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Score") +
  ylim(c(-20, 30))+
  ggtitle("Gamma Score")

##---
alpha_l <- lapply(b_est$model.stats$optimizer$updates$alpha[[1]][58:79], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
alpha <- do.call(rbind, Map(cbind, it = 57 + seq_along(alpha_l), alpha_l))
p_score_a <- ggplot(data = alpha,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Score") +
  ggtitle("Alpha Score")


##---
mufpc1_l <- lapply(b_est$model.stats$optimizer$updates$mu[[1]], 
                  function(x) {
                    data.frame(
                      "par" = paste0("par", seq_along(x$xscore)),
                      "s_like" = x$xscore0, 
                      "s_prio" = x$xscore - x$xscore0, 
                      "s_post" = x$xscore) 
                  })
mufpc1 <- do.call(rbind, Map(cbind, it = seq_along(mufpc1_l), mufpc1_l))
p1 <- ggplot(data = mufpc1,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-10, 10))+
  ggtitle("Mu Fpc1 Likelihood Score")
p2 <- ggplot(data = mufpc1,
       aes(x = it, y = s_prio, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Fpc1 Priori Score")
grid.arrange(p1, p2, nrow = 1)
p_score_fpc1 <- ggplot(data = mufpc1,
             aes(x = it, y = s_post, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Score") +
  ylim(c(-10, 10)) +
  ggtitle("Mu Fpc1 Score")
grid.arrange(p1, p2, nrow = 1)

mufpc2_l <- lapply(b_est$model.stats$optimizer$updates$mu[[2]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
mufpc2 <- do.call(rbind, Map(cbind, it = seq_along(mufpc2_l), mufpc2_l))
ggplot(data = mufpc2,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-10, 10))+
  ggtitle("Mu Fpc2 Score")

mufpc3_l <- lapply(b_est$model.stats$optimizer$updates$mu[[3]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
mufpc3 <- do.call(rbind, Map(cbind, it = seq_along(mufpc3_l), mufpc3_l))
ggplot(data = mufpc3,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-10, 10))+
  ggtitle("Mu Fpc3 Score")

mufpc4_l <- lapply(b_est$model.stats$optimizer$updates$mu[[4]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
mufpc4 <- do.call(rbind, Map(cbind, it = seq_along(mufpc4_l), mufpc4_l))
ggplot(data = mufpc4,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-10, 10))+
  ggtitle("Mu Fpc4 Score")

mufpc5_l <- lapply(b_est$model.stats$optimizer$updates$mu[[5]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
mufpc5 <- do.call(rbind, Map(cbind, it = seq_along(mufpc5_l), mufpc5_l))
ggplot(data = mufpc5,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-10, 10))+
  ggtitle("Mu Fpc5 Score")

mup_l <- lapply(b_est$model.stats$optimizer$updates$mu[[6]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "s_like" = x$xscore0, 
                       "s_prio" = x$xscore - x$xscore0, 
                       "s_post" = x$xscore) 
                   })
mup <- do.call(rbind, Map(cbind, it = seq_along(mup_l), mup_l))
ggplot(data = mup,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-5000, 7500))+
  ggtitle("Mu P Score")
p_score_mup <- ggplot(data = mup,
       aes(x = it, y = s_post, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Score") +
  ylim(c(-5000, 7500))+
  ggtitle("Mu P Score")


##---
sigma_l <- lapply(b_est$model.stats$optimizer$updates$sigma[[1]], 
                function(x) {
                  data.frame(
                    "par" = paste0("par", seq_along(x$xscore)),
                    "s_like" = x$xscore0, 
                    "s_prio" = x$xscore - x$xscore0, 
                    "s_post" = x$xscore) 
                })
sigma <- do.call(rbind, Map(cbind, it = seq_along(sigma_l), sigma_l))
ggplot(data = sigma,
       aes(x = it, y = s_like, col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ylim(c(-2500, 8000)) +
  ggtitle("Sigma Score")



# Hesse Trace Plots -------------------------------------------------------

lambda_l <- lapply(b_est$model.stats$optimizer$updates$lambda[[2]], 
                   function(x) {
                     data.frame(
                       "h_like" = det(x$xh0)^-1, 
                       "h_prio" = det(x$xh - x$xh0)^-1, 
                       "h_post" = det(x$xh)^-1)
                   })
lambda <- do.call(rbind, Map(cbind, it = seq_along(lambda_l), lambda_l))
ggplot(data = lambda,
       aes(x = it, y = log(h_like))) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "log(|-H|^-1)") +
 # ylim(c(0, 10000000)) +
  ggtitle("Lambda Inv(Det(Neg(Hess)))")

p_hess_l <- ggplot(data = lambda,
       aes(x = it, y = h_post)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  # ylim(c(0, 10000000)) +
  ggtitle("Lambda Inv(Det(Neg(Hess)))")
grid.arrange(p_score_l, p_hess_l, nrow = 1)

##---
gamma_l <- lapply(b_est$model.stats$optimizer$updates$gamma[[1]], 
                  function(x) {
                    data.frame(
                      "h_like" = det(x$xh0)^-1, 
                      "h_prio" = det(x$xh - x$xh0)^-1, 
                      "h_post" = det(x$xh)^-1)
                  })
gamma <- do.call(rbind, Map(cbind, it = seq_along(gamma_l), gamma_l))
ggplot(data = gamma,
       aes(x = it, y = h_like)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  ggtitle("Gamma Inv(Det(Neg(Hess)))")

p_hess_g <- ggplot(data = gamma,
       aes(x = it, y = h_post)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  ggtitle("Gamma Inv(Det(Neg(Hess)))")
grid.arrange(p_score_g, p_hess_g, nrow = 1)

##---
alpha_l <- lapply(b_est$model.stats$optimizer$updates$alpha[[1]][58:80], 
                  function(x) {
                    data.frame(
                      "h_like" = det(x$xh0)^-1, 
                      "h_prio" = det(x$xh - x$xh0)^-1, 
                      "h_post" = det(x$xh)^-1)
                  })
alpha <- do.call(rbind, Map(cbind, it = 57 + seq_along(alpha_l), alpha_l))
ggplot(data = alpha,
       aes(x = it, y = h_like)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  ggtitle("Alpha Inv(Det(Neg(Hess)))")
p_hess_a <- ggplot(data = alpha,
       aes(x = it, y = h_post)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  ggtitle("Alpha Inv(Det(Neg(Hess)))")
grid.arrange(p_score_a, p_hess_a, nrow = 1)



##---
mufpc1_l <- lapply(b_est$model.stats$optimizer$updates$mu[[1]], 
                   function(x) {
                     data.frame(
                       "par" = paste0("par", seq_along(x$xscore)),
                       "h_like" = diag(x$xh0)^-1, 
                       "h_prio" = diag(x$xh - x$xh0)^-1, 
                       "h_post" = diag(x$xh)^-1)
                   })
mufpc1 <- do.call(rbind, Map(cbind, it = seq_along(mufpc1_l), mufpc1_l))
ggplot(data = mufpc1,
       aes(x = it, y = log(h_like), col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "log(diag(-H)^-1)") +
  ggtitle("Mu FPC1 Inv(Diag(Neg(Hess)))")
p_hess_fpc1 <- ggplot(data = mufpc1,
       aes(x = it, y = log(h_post), col = par)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "log(diag(-H)^-1)") +
  ggtitle("Mu FPC1 Inv(Diag(Neg(Hess)))")
grid.arrange(p_score_a, p_hess_a, nrow = 1)

mup_l <- lapply(b_est$model.stats$optimizer$updates$mu[[6]], 
                  function(x) {
                    data.frame(
                      "h_like" = det(x$xh0)^-1, 
                      "h_prio" = det(x$xh - x$xh0)^-1, 
                      "h_post" = det(x$xh)^-1)
                  })
mup <- do.call(rbind, Map(cbind, it = seq_along(mup_l), mup_l))
p_hess_mup <- ggplot(data = mup,
                   aes(x = it, y = h_post)) +
  geom_line() +
  geom_vline(xintercept = 57, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "|-H|^-1") +
  ggtitle("Mu P Inv(Det(Neg(Hess)))")
grid.arrange(p_score_mup, p_hess_mup, nrow = 1)



# Remove Subs with Only 1 Obs ---------------------------------------------

n1obs <- which(table(d_rirs$marker, d_rirs$i)[1, ] != 1)

f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_along(4), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, 
                data = d_rirs_tru %>% filter(id %in% n1obs), 
                timevar = "year", maxit = 1500, n.iter = 900,
                burnin = 1000, thin = 3, verbose = TRUE)





# RI+RS Model -------------------------------------------------------------


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
                burnin = 1000, thin = 3, verbose = FALSE)




# Differences Simulation --------------------------------------------------

# Covariance matrix own simulation
auto <- matrix(c(0.08, -0.07, -0.07, 0.9), ncol = 2)
cross <- matrix(rep(0.03, 4), ncol = 2)
cor <- matrix(c(0, 1, 0.75, 0.5, 0, 0,
                1, 0, 1, 0.75, 0.5, 0,
                0.75, 1, 0, 1, 0.75, 0.5,
                0.5, 0.75, 1, 0, 1, 0.75,
                0, 0.5, 0.75, 1, 0, 1,
                0, 0, 0.5, 0.75, 1, 0),
              ncol = 6)
cov <- kronecker(cor, cross) + kronecker(diag(c(1, 1.2, 1.4, 1.6, 1.8, 2)), 
                                         auto)

upper <- round(cov, digits = 2)
upper_cor <- round(cov2cor(cov), digits = 2)

# For printing in xtable
upper[upper.tri(cov)] <- ""
upper_cor[upper.tri(cov)] <- ""
upper <- as.data.frame(upper)
upper_cor <- as.data.frame(upper_cor)
print(xtable(upper), include.rownames = FALSE, include.colnames = FALSE)
print(xtable(upper_cor), include.rownames = FALSE, include.colnames = FALSE)


# Mauf Simulation
D_on1 <- D
D_on1[seq(2, 12, by = 2), ] <- D_on1[seq(2, 12, by = 2), ] * 25
D_on1[, seq(2, 12, by = 2)] <- D_on1[, seq(2, 12, by = 2)] * 25
upper <- as.matrix(round(D_on1, digits = 2))
upper_cor <- as.matrix(round(cov2cor(D_on1), digits = 2))

# For printing in xtable
upper[upper.tri(cov)] <- ""
upper_cor[upper.tri(cov)] <- ""
upper <- as.data.frame(upper)
upper_cor <- as.data.frame(upper_cor)
print(xtable(upper), include.rownames = FALSE, include.colnames = FALSE)
print(xtable(upper_cor), include.rownames = FALSE, include.colnames = FALSE)



# Multivariate Failure Model ----------------------------------------------

# Model is from scenarioI_mauff_nu.R on sunflower
# i <- 2
# nfpc <- 12
# b_est <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru,
#                 timevar = "year", maxit = 1480, par_trace = TRUE,
#                 sampler = FALSE, update_nu = TRUE, verbose = TRUE)
b_est <- readRDS(paste0(server_wd, "scen_mauff/bamlss_tru_nu/b_2_full_fail.Rds"))
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:1480) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_vline(xintercept = 31, linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters")
