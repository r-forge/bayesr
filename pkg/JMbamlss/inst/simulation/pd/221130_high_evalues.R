location <- "workstation"


if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


if (location == "server_linux") {
  results_wd <- paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
                       "simulation/")
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
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
devtools::load_all()


# Calculate the MFPCs on each simulation run ------------------

data_fpc_i <- function(i, dat_setting = "scen_I_130922", ...) {
  
  # Load the data
  load(paste0(results_wd, dat_setting, "/data/d", i, ".Rdata"))
  
  # Estimate the model using estimated FPCs
  preproc_MFPCA(d_rirs$data,
                uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                npc = 2, nbasis = 4, ...)
  
}
data_fpc <- Vectorize(data_fpc_i, vectorize.args = "i", SIMPLIFY = FALSE)
mfpc_est2 <- data_fpc(100:199, remove_obs = 2, save_uniFPCA = TRUE)



# Find simulation iterations with high eigenvalue -------------------------

# Simulation runs with high multivariate eigenvalues
evals2 <- sapply(mfpc_est2, "[[", "values")
high_ev2 <- which(evals2[1, ] > 10)

# Extract univariate eigenvalues for these iterations
all_ev <- sapply(mfpc_est2, function (i) {
  sapply(attr(i, "uniFPC"), "[[", "evalues")
})
all_ev[, high_ev2]
summary(c(all_ev[, -high_ev2]))
summary(c(all_ev[, high_ev2]))
summary(colSums(all_ev[, -high_ev2]))
summary(colSums(all_ev[, high_ev2]))

# Maximal multivariate scores of iteration runs with high multivariate evalues
lapply(high_ev2, function (i) {
  max(abs(mfpc_est2[[i]]$scores[, 1]))
})

# Maximal scores over all the FPCs
lapply(high_ev2, function (i) {
  apply(mfpc_est2[[i]]$scores, 2, function(m) max(abs(m)))
})

# Number of longitudinal observations for high-score observations
cut_off <- 10
high_sc <- lapply(mfpc_est2, function(it) {
  obs <- which(apply(it$scores, 1, function (sc) any(abs(sc) > cut_off)))
  if (length(obs)) {
    FPCA <- attr(it, "uniFPCA")
    # Sollte man für jedes obs getrennt machen!
    do.call(rbind, lapply(obs, function(id) {
      sc_obs <- sapply(FPCA, function (m) {
        m$scores[id, ]
      })
      m_obs <- which(apply(sc_obs, 2, function(v) any(abs(v) > cut_off)))
      if (length(m_obs)) {
        m_fpc <- sapply(m_obs, function (m) {
          paste0("fpc", which(abs(FPCA[[m]]$scores[id, ]) > cut_off))
        })
        m_sco <- sapply(m_obs, function (m) {
          FPCA[[m]]$scores[id, which(abs(FPCA[[m]]$scores[id, ]) > cut_off)]
        })
        n_obs <- sapply(m_obs, function (m) {
          sum(!is.na(FPCA[[m]]$Y[id, ]))
        })
        if (is.list(m_fpc)) {
          n_obs <- rep(n_obs, sapply(m_fpc, length))
          m_fpc <- unlist(m_fpc)
          m_sco <- unlist(m_sco)
        }
        data.frame(id = id, marker = names(n_obs), n_obs = as.vector(n_obs), 
                   fpc = as.vector(m_fpc), score = as.vector(m_sco))
      } else {
        data.frame(id = id, marker = NA, n_obs = NA, 
                   fpc = NA, score = NA)
      }
    }))
  } else {
    data.frame(id = NA, marker = NA, n_obs = NA, 
               fpc = NA, score = NA)
  }
})
high_sc <- do.call(rbind, Map(cbind, it = seq_along(high_sc), high_sc))
high_sc <- high_sc[!is.na(high_sc$id), ]

table(high_sc$fpc)
ggplot(na.omit(high_sc), aes(x = as.factor(n_obs), y = log(abs(score)))) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Number of Observations for Multivariate Scores > 10", 
          "Most multivariate scores come from univariate scores > 10") +
  labs(x = "Number of Longitudinal Observations",
       y = "Log(abs(Univariate Score))")



# Look at Specific Iteration with 2 high Scores ---------------------------

# Pick Iteration 111 because the first two FPCs have high multivariate scores
m <- mfpc_est2[[high_ev2[[3]]]]
m$scores[, 1:2]

FPCA <- attr(m, "uniFPCA")
w <- which.min(FPCA$m2$scores[, 1])

dat <- as.data.frame(FPCA$m2$Y) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = 1:100, names_to = "time") %>%
  na.omit() %>%
  mutate(time = as.numeric(time),
         idcol = ifelse(row_id == w, TRUE, FALSE)) %>%
  arrange(idcol)
ggplot(dat %>% filter(idcol == FALSE),
       aes(x = time, y = value, group = row_id, color = idcol)) +
  geom_line(alpha = 0.7) +
  geom_line(data = dat %>% filter(idcol == TRUE)) +
  theme_bw() +
  ggtitle("Simulation 111: High Univariate Score", 
          "Evalue1: 0.379, Score41: -654") +
  labs(x = "Time", y = "Centered Obs", color = element_blank()) +
  scale_color_manual(labels = c("Rest", "ID41"), values = c("grey", "red"))


w <- which.max(FPCA$m5$scores)
dat <- as.data.frame(FPCA$m5$Y) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = 1:98, names_to = "time") %>%
  na.omit() %>%
  mutate(time = as.numeric(time),
         idcol = ifelse(row_id == w, TRUE, FALSE)) %>%
  arrange(idcol)
ggplot(dat %>% filter(idcol == FALSE),
       aes(x = time, y = value, group = row_id, color = idcol)) +
  geom_line(alpha = 0.7) +
  geom_line(data = dat %>% filter(idcol == TRUE)) +
  theme_bw() +
  ggtitle("Simulation 111: High Univariate Score", 
          "Evalue2: 0.379, Score27: 32.33") +
  labs(x = "Time", y = "Centered Obs", color = element_blank()) +
  scale_color_manual(labels = c("Rest", "ID27"), values = c("grey", "red"))


# Remove all observations with less than 3 observations -------------------


mfpc_est3 <- data_fpc(100:199, remove_obs = 3, save_uniFPCA = TRUE)
evals3 <- sapply(mfpc_est3, "[[", "values")
high_ev3 <- which(evals3[1, ] > 10)

m <- mfpc_est3[[high_ev3]]
w <- which.min(m$scores[, 1])
FPCA <- attr(m, "uniFPCA")
mFData <- multiFunData(lapply(FPCA, function (mark) {
  funData(argvals = mark$argvals, X = mark$Yhat)
}))
uni_scores_w <- sapply(FPCA, function (x){
  x$scores[w, ]
})

# Plot the observations
dat <- as.data.frame(FPCA$m5$Y) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -ncol(.), names_to = "time") %>%
  na.omit() %>%
  mutate(time = as.numeric(time),
         idcol = ifelse(row_id == w, TRUE, FALSE)) %>%
  arrange(idcol)
ggplot(dat %>% filter(idcol == FALSE),
       aes(x = time, y = value, group = row_id, color = idcol)) +
  geom_line(alpha = 0.7) +
  geom_line(data = dat %>% filter(idcol == TRUE)) +
  theme_bw() +
  ggtitle("Simulation 155: High Univariate Score (After Removing #Obs<3)", 
          "Evalue1: 0.484, Score82: -35.90") +
  labs(x = "Time", y = "Centered Obs", color = element_blank()) +
  scale_color_manual(labels = c("Rest", "ID82"), values = c("grey", "red"))

# Plot the fitted curve
dat <- as.data.frame(FPCA$m5$Yhat) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -ncol(.), names_to = "time") %>%
  na.omit() %>%
  mutate(time = as.numeric(time),
         idcol = ifelse(row_id == w, TRUE, FALSE)) %>%
  arrange(idcol)
ggplot(dat %>% filter(idcol == FALSE),
       aes(x = time, y = value, group = row_id, color = idcol)) +
  geom_line(alpha = 0.7) +
  geom_line(data = dat %>% filter(idcol == TRUE)) +
  theme_bw() +
  ggtitle("Simulation 155: High Univariate Score (After Removing #Obs<3)", 
          "Evalue1: 0.484, Score82: -35.90") +
  labs(x = "Time", y = "Centered Obs", color = element_blank()) +
  scale_color_manual(labels = c("Rest", "ID82"), values = c("grey", "red"))
