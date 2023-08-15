
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
                                           "simulation"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation",
                    "server_windows" = "H:/JMbamlss/simulation")


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



# Extract FPCs from models ------------------------------------------------


extract_fpc_info_i <- function (i, results_wd, setting, folder) {
  
  b_est <- readRDS(file.path(results_wd, setting, folder, i))
  
  list("nfpc" = attr(b_est, "nfpc"),
       "mfpca" = attr(b_est, "FPCs"))
}
extract_fpc_info <- Vectorize(extract_fpc_info_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
fpc_info <- function (results_wd, setting, folder) {
  
  m <- list.files(file.path(results_wd, setting, folder))
  extract_fpc_info(m, results_wd = results_wd, setting = setting,
                   folder = folder)
  
}


simI <- fpc_info(results_wd = server_wd,
                 setting = "scen_I_230719",
                 folder = "bamlss_est95")
simII <- fpc_info(results_wd = server_wd,
                  setting = "scen_II_230719",
                  folder = "bamlss_est95")



# Number of FPCs in the Model ---------------------------------------------

summary(sapply(simI, "[[", "nfpc"))
summary(sapply(simII, "[[", "nfpc"))


# Extract Eigenfunctions and Combine to Data.frame ------------------------

simI_tru <- extract_fpc_info_i("b100.rds", server_wd, "scen_I_230719", 
                               "bamlss_tru")
simII_tru <- extract_fpc_info_i("b100.rds", server_wd, "scen_II_230719", 
                                "bamlss_tru")



# Calculate the norms for evaluation of estimation accuracy ---------------


eigenfun_norm_i <- function(mfpca, mfpca_tru) {
  
  # Check whether the argvals are the same
  if (!isTRUE(all.equal(argvals(mfpca$functions), 
                        argvals(mfpca_tru$functions)))) {
    argvals <- argvals(mfpca_tru$functions)[[1]][[1]]
    new_mfpcs <- JMbamlss:::eval_mfpc(mfpca, argvals)
    new_mfpcs <- unname(
      lapply(split(new_mfpcs, as.factor(rownames(new_mfpcs))), function(x) {
        funData(argvals, X = matrix(x, ncol = length(argvals), byrow = TRUE))
    }))
    mfpca$functions <- multiFunData(new_mfpcs)
  }
  
  sapply(seq_len(nObs(mfpca$functions)), function (fu_it) {
    tr <- extractObs(mfpca_tru$functions, fu_it)
    es <- flipFuns(tr, extractObs(mfpca$functions, fu_it))
    norm(tr - es)
  })
}

eigenfun_norm <- Vectorize(eigenfun_norm_i, vectorize.args = "mfpca",
                           SIMPLIFY = TRUE)



norm_simI <- eigenfun_norm(lapply(simI, "[[", "mfpca"), simI_tru$mfpca)
norm_simII <- eigenfun_norm(lapply(simII, "[[", "mfpca"), simII_tru$mfpca)


# Plot the estimated MFPCs ------------------------------------------------


mfpc2dataframe_i <- function (mfpca, mfpca_tru) {
  
  # Check whether the argvals are the same
  if (!isTRUE(all.equal(argvals(mfpca$functions), 
                        argvals(mfpca_tru$functions)))) {
    argvals <- argvals(mfpca_tru$functions)[[1]][[1]]
    new_mfpcs <- JMbamlss:::eval_mfpc(mfpca, argvals)
    new_mfpcs <- unname(
      lapply(split(new_mfpcs, as.factor(rownames(new_mfpcs))), function(x) {
        funData(argvals, X = matrix(x, ncol = length(argvals), byrow = TRUE))
      }))
    mfpca$functions <- multiFunData(new_mfpcs)
  }
  
  mfpca$functions <- flipFuns(mfpca_tru$functions, mfpca$functions)
  
  dat <- lapply(mfpca$functions@.Data, function(m) {
    data.frame(t = m@argvals[[1]], 
               X = c(t(m@X)), 
               fpc = factor(paste0("fpc.", rep(seq_along(mfpca$values), 
                                               each = length(m@argvals[[1]]))),
                            levels = paste0("fpc.", seq_along(mfpca$values))))
  })
  do.call(rbind, Map(cbind, marker = factor(paste0("m", seq_along(dat))), dat))
}


mfpc2dataframe <- Vectorize(mfpc2dataframe_i, SIMPLIFY = FALSE, 
                            vectorize.args = "mfpca")

estimated_mfpcs <- function(mfpca_list, mfpca_tru) {
  dat <- mfpc2dataframe(mfpca_list, mfpca_tru)
  do.call(rbind, Map(cbind, it = factor(seq_along(dat)), dat))
}

# Simulation Scenario I
dat_mfpc_est <- estimated_mfpcs(lapply(simI, "[[", "mfpca"), 
                                mfpca_tru = simI_tru$mfpca)
dat_mfpc_tru <- mfpc2dataframe_i(mfpca = simI_tru$mfpca,
                                 mfpca_tru = simI_tru$mfpca)

ggplot(dat_mfpc_est %>% filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3")) %>%
         mutate(marker = factor(marker, labels = paste("Outcome", 1:6)),
                fpc = factor(droplevels(fpc), labels = paste("MFPC", 1:3))),
       aes(x = t, y = X, group = it)) +
  geom_line(alpha = 0.1) +
  geom_line(data = dat_mfpc_tru %>% 
              filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3")) %>%
              mutate(it = 500,
                     marker = factor(marker, labels = paste("Outcome", 1:6)),
                     fpc = factor(droplevels(fpc), 
                                  labels = paste("MFPC", 1:3))),
            color = "red") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1), labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(fpc ~ marker) +
  labs(y = NULL, x = "Time")
# Save 6x12

# Simulation Scenario II
dat_mfpc_est <- estimated_mfpcs(lapply(simII, "[[", "mfpca"), 
                                mfpca_tru = simII_tru$mfpca)
dat_mfpc_tru <- mfpc2dataframe_i(mfpca = simII_tru$mfpca,
                                 mfpca_tru = simII_tru$mfpca)

p1 <- ggplot(dat_mfpc_est %>% filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3")) %>%
         mutate(marker = factor(marker, labels = paste("Outcome", 1:2)),
                fpc = factor(droplevels(fpc), labels = paste("MFPC", 1:3))),
       aes(x = t, y = X, group = it)) +
  geom_line(alpha = 0.1) +
  geom_line(data = dat_mfpc_tru %>% 
              filter(fpc %in% c("fpc.1", "fpc.2", "fpc.3")) %>%
              mutate(it = 500,
                     marker = factor(marker, labels = paste("Outcome", 1:2)),
                     fpc = factor(droplevels(fpc), 
                                  labels = paste("MFPC", 1:3))),
            color = "red") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1), labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(fpc ~ marker) +
  labs(y = NULL, x = "Time")
p2 <- ggplot(dat_mfpc_est %>% filter(fpc %in% c("fpc.4", "fpc.5", "fpc.6")) %>%
         mutate(marker = factor(marker, labels = paste("Outcome", 1:2)),
                fpc = factor(droplevels(fpc), labels = paste("MFPC", 4:6))),
       aes(x = t, y = X, group = it)) +
  geom_line(alpha = 0.1) +
  geom_line(data = dat_mfpc_tru %>% 
              filter(fpc %in% c("fpc.4", "fpc.5", "fpc.6")) %>%
              mutate(it = 500,
                     marker = factor(marker, labels = paste("Outcome", 1:2)),
                     fpc = factor(droplevels(fpc), 
                                  labels = paste("MFPC", 4:6))),
            color = "red") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1), labels = c(0, 0.5, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(fpc ~ marker) +
  labs(y = NULL, x = "Time")
gridExtra::grid.arrange(p1, p2, nrow = 1)
# Save 4x10