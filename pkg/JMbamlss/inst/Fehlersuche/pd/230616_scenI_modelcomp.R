library(tidyverse)

# Convenience function
acc <- function(model){
  summ <- summary(model$samples[[1]][, grep("accepted", 
                                            colnames(model$samples[[1]]))])
  matrix(summ$statistics[, 1], ncol = 1, 
         dimnames = list(names(summ$statistics[, 1]), "Mean"))
}



server_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.",
                    "hu-berlin.de,share=volkmana.hub/JMbamlss/simulation")



# Model After Code Updates ------------------------------------------------


load(file.path(server_wd, "scen_I_230615/bamlss_est", "b100.Rdata"))
acc(b_est)
xtable::xtable(acc(b_est))

# Variance estimates for leading FPCs
fpc1 <- grep("\\(id,fpc\\.1\\)", colnames(b_est$samples[[1]]))[1:150]
fpc2 <- grep("\\(id,fpc\\.2\\)", colnames(b_est$samples[[1]]))[1:150]
all_draws <- data.frame(
  id = rep(rep(seq_len(150), each = 1001), times = 2),
  re = c(b_est$samples[[1]][, fpc1], b_est$samples[[1]][, fpc2]),
  type = factor(rep(c("FPC 1", "FPC 2"), each = 1001*150))
)
fpc1_tau2 <- unname(summary(b_est$samples[[1]][, fpc1[150] + 1])$statistics[1])
fpc2_tau2 <- unname(summary(b_est$samples[[1]][, fpc2[150] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-3, 3, by = 0.1), seq(-3, 3, by = 0.1)),
  y = c(dnorm(seq(-3, 3, by = 0.1), 0, sqrt(fpc1_tau2)),
        dnorm(seq(-3, 3, by = 0.1), 0, sqrt(fpc2_tau2))),
  type = factor(c(rep("FPC 1", length(seq(-3, 3, by = 0.1))),
                  rep("FPC 2", length(seq(-3, 3, by = 0.1)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("NEW MODELS -- Distribution of All FPC Score Samples", 
          paste0("MJM ScenI (bamlss), 150 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Scores", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")

# Variance estimates for trailing FPCs
fpc6 <- grep("\\(id,fpc\\.6\\)", colnames(b_est$samples[[1]]))[1:150]
fpc7 <- grep("\\(id,fpc\\.7\\)", colnames(b_est$samples[[1]]))[1:150]
all_draws <- data.frame(
  id = rep(rep(seq_len(150), each = 1001), times = 2),
  re = c(b_est$samples[[1]][, fpc6], b_est$samples[[1]][, fpc7]),
  type = factor(rep(c("FPC 6", "FPC 7"), each = 1001*150))
)
fpc6_tau2 <- unname(summary(b_est$samples[[1]][, fpc6[150] + 1])$statistics[1])
fpc7_tau2 <- unname(summary(b_est$samples[[1]][, fpc7[150] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-2.5, 2.5, by = 0.1), seq(-1, 1, by = 0.05)),
  y = c(dnorm(seq(-2.5, 2.5, by = 0.1), 0, sqrt(fpc6_tau2)),
        dnorm(seq(-1, 1, by = 0.05), 0, sqrt(fpc7_tau2))),
  type = factor(c(rep("FPC 6", length(seq(-2.5, 2.5, by = 0.1))),
                  rep("FPC 7", length(seq(-1, 1, by = 0.05)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("NEW MODELS -- Distribution of All FPC Score Samples", 
          paste0("MJM ScenI (bamlss), 150 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Scores", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")

# Trace plot of updating the alpha parameter
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:418) %>%
         pivot_longer(cols = -it) %>%
         filter(it %in% 142:418),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(1.5, 0.6, 0.3, -0.3, -0.6, -1.5),
             linetype = "dotted") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Model 100: Updating Trace of Alpha Parameters")
rm(b_est)


# Model Before Code Update ------------------------------------------------


load(file.path(server_wd, "scen_I_221019/bamlss_est", "b100.Rdata"))
acc(b_est)
# xtable::xtable(acc(b_est))

# Variance estimates for leading FPCs
fpc1 <- grep("\\(id,fpc\\.1\\)", colnames(b_est$samples[[1]]))[1:150]
fpc2 <- grep("\\(id,fpc\\.2\\)", colnames(b_est$samples[[1]]))[1:150]
all_draws <- data.frame(
  id = rep(rep(seq_len(150), each = 1001), times = 2),
  re = c(b_est$samples[[1]][, fpc1], b_est$samples[[1]][, fpc2]),
  type = factor(rep(c("FPC 1", "FPC 2"), each = 1001*150))
)
fpc1_tau2 <- unname(summary(b_est$samples[[1]][, fpc1[150] + 1])$statistics[1])
fpc2_tau2 <- unname(summary(b_est$samples[[1]][, fpc2[150] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-3, 3, by = 0.1), seq(-3, 3, by = 0.1)),
  y = c(dnorm(seq(-3, 3, by = 0.1), 0, sqrt(fpc1_tau2)),
        dnorm(seq(-3, 3, by = 0.1), 0, sqrt(fpc2_tau2))),
  type = factor(c(rep("FPC 1", length(seq(-3, 3, by = 0.1))),
                  rep("FPC 2", length(seq(-3, 3, by = 0.1)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("OLD MODELS -- Distribution of All FPC Score Samples", 
          paste0("MJM ScenI (bamlss), 150 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Scores", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")

# Variance estimates for trailing FPCs
fpc6 <- grep("\\(id,fpc\\.6\\)", colnames(b_est$samples[[1]]))[1:150]
fpc7 <- grep("\\(id,fpc\\.7\\)", colnames(b_est$samples[[1]]))[1:150]
all_draws <- data.frame(
  id = rep(rep(seq_len(150), each = 1001), times = 2),
  re = c(b_est$samples[[1]][, fpc6], b_est$samples[[1]][, fpc7]),
  type = factor(rep(c("FPC 6", "FPC 7"), each = 1001*150))
)
fpc6_tau2 <- unname(summary(b_est$samples[[1]][, fpc6[150] + 1])$statistics[1])
fpc7_tau2 <- unname(summary(b_est$samples[[1]][, fpc7[150] + 1])$statistics[1])
all_taus <- data.frame(
  x = c(seq(-2.5, 2.5, by = 0.1), seq(-1, 1, by = 0.05)),
  y = c(dnorm(seq(-2.5, 2.5, by = 0.1), 0, sqrt(fpc6_tau2)),
        dnorm(seq(-1, 1, by = 0.05), 0, sqrt(fpc7_tau2))),
  type = factor(c(rep("FPC 6", length(seq(-2.5, 2.5, by = 0.1))),
                  rep("FPC 7", length(seq(-1, 1, by = 0.05)))))
)
ggplot(all_draws, aes(x = re)) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~type, scales = "free") +
  ggtitle("OLD MODELS -- Distribution of All FPC Score Samples", 
          paste0("MJM ScenI (bamlss), 150 Ids, 1001 Samples, ",
                 "Normal with Posterior Mean of tau2")) +
  labs(x = "Random Scores", y = "Density") +
  geom_line(data = all_taus, aes(x = x, y = y), linetype = "dashed")
rm(b_est)


# Model Iteration with Previous Numerical Problems  -----------------------

load(file.path(server_wd, "scen_I_230615/bamlss_est", "b151.Rdata"))
acc(b_est)

# Trace plot of updating the alpha parameter
alpha <- sapply(b_est$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:498) %>%
         pivot_longer(cols = -it) %>%
         filter(it %in% 400:498),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = c(1.5, 0.6, 0.3, -0.3, -0.6, -1.5),
             linetype = "dotted") +
  theme_bw() +
  labs(x = "Iteration", y =  "Estimate", col = NULL) +
  ggtitle("Model 151: Updating Trace of Alpha Parameters")
rm(b_est)
