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
source("R/sim_helperfun.R")
source("R/preprocessing.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/fpca.R")


load(paste0(results_wd, "scen_I_130922/bamlss_est/b102.Rdata"))
names(attributes(b_est))
MFPCs <- attr(b_est, "FPCs")


head(colnames(samples(b_est, model = "mu")), n = 155)
summary(samples(b_est)[ ,grep("mu\\.s\\.s\\(id,fpc\\.[0-9]+\\)\\.tau21", 
                              colnames(samples(b_est)))])

# Wie viel Prozent der Varianz in den Daten erklären die FPCs?
# Was sind die geschätzten Eigenwerte?
# Sind die gezogenen Random Effects unkorreliert?


# Extract samples and number of FPCs of models -----------------
extract_fpc_info_i <- function (i, setting, folder) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  
  list("nfpc" = attr(b_est, "nfpc"),
       "samples" = samples(b_est, model = "mu"))
}
extract_fpc_info <- Vectorize(extract_fpc_info_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
fpc_info <- function (setting, folder) {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_fpc_info(m, setting = setting, folder = folder)
  
}

b_est_1 <- fpc_info(setting = "scen_I_221019", folder = "bamlss_est_1")
b_est_975 <- fpc_info(setting = "scen_I_221019", folder = "bamlss_est")


# Calculate the MFPCs on each simulation run ------------------

data_fpc_i <- function(i, dat_setting = "scen_I_130922") {

   # Load the data
  load(paste0(results_wd, dat_setting, "/data/d", i, ".Rdata"))
  
  # Estimate the model using estimated FPCs
  few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1, 
                   function (x) any(x < 2))
  mfpca_est <- preproc_MFPCA(d_rirs$data %>%
                               filter(id %in% paste(which(!few_obs))) %>% 
                               droplevels(), 
                             uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                             npc = 2, nbasis = 4)

  mfpca_est
}
data_fpc <- Vectorize(data_fpc_i, vectorize.args = "i", SIMPLIFY = FALSE)
mfpc_est <- data_fpc(100:199)

# For comparison use the true calculated MFPC
load(paste0(results_wd, "scen_I_130922/bamlss_tru/b102.Rdata"))

# Compare Eigenvalues -----------------------------------------
eval_dat <- data.frame(do.call(rbind, lapply(mfpc_est, "[[", "values"))) %>%
  rbind(attr(b_est, "FPCs")$values) %>%
  pivot_longer(X1:X12, names_to = "Name", values_to = "Eigenvalue") %>%
  mutate(Name = factor(Name, levels = paste0("X", 1:12), 
                       labels = paste0("FPC", 1:12)),
         src = factor(rep(c("Estimated", "True"), times = c(100*12, 12))),
         it = rep(1:101, each = 12))
ggplot(eval_dat, aes(x = Name, y = log(Eigenvalue), color = src)) +
  geom_boxplot() +
  scale_color_manual(values = c("black", "red"), name = NULL) +
  ggtitle("Estimated/True Eigenvalues in the Simulation Scenario I")


# Compare univariate norms of the MFPCs ------------------------------
norms_sq <- lapply(mfpc_est, function (x){
  sapply(x$functions@.Data, funData::norm, squared = TRUE)
})

norms_dat <- do.call(rbind, norms_sq) %>%
  rbind(sapply(attr(b_est,"FPCs")$functions@.Data, funData::norm,
               squared = TRUE)) %>%
  data.frame() %>%
  mutate(FPC = factor(rep(paste0("FPC", 1:12), times = 101), 
                      levels = paste0("FPC", 1:12)),
         src = factor(rep(c("Estimated", "True"), times = c(100*12, 12))),
         it = rep(1:101, each = 12)) %>%
  pivot_longer(X1:X6, names_to = "Marker", values_to = "UniNorms") %>%
  mutate(Marker = factor(Marker, levels = paste0("X", 1:6), 
                         labels = paste0("m", 1:6)))

ggplot(norms_dat, aes(x = Marker, y = UniNorms, group = it, color = src)) +
  theme_bw() +
  geom_line(aes(alpha = src)) +
  scale_color_manual(values = c("black", "red"), name = NULL) +
  scale_alpha_manual(values = c(0.1, 1)) +
  guides(alpha = "none") +
  facet_wrap(vars(FPC), nrow = 4) +
  ggtitle("Univariate Norms of FPCs in the Simulation Scenario I")


# Are the MFPCs of the problematic simulated data different? -------------

m_95 <- list.files(paste0(results_wd, "scen_I_130922/bamlss_est/"))
m_975 <- list.files(paste0(results_wd, "scen_I_221019/bamlss_est/"))

which(m_95 == "b170.Rdata")

match(m_95[seq_len(which(m_95 == "b170.Rdata"))], m_975)
match(m_975, m_95[seq_len(which(m_95 == "b170.Rdata"))])


# Are the preliminary MFPCs different for the problematic data -------------

probl <- which(! paste0("b", 100:199, ".Rdata") %in% m_95)

ggplot(eval_dat %>% filter(src == "Estimated"), 
       aes(x = Name, y = log(Eigenvalue))) +
  geom_boxplot() +
  theme_bw() +
  geom_point(data = eval_dat %>% filter(it %in% probl),
             aes(x = Name, y = log(Eigenvalue)), shape = 8, color = "orange") +
  ggtitle("Estimated Eigenvalues for Problematic Iterations")

ggplot(norms_dat %>% filter(src == "Estimated") %>%
         mutate(problem = factor(it %in% probl)),
       aes(x = Marker, y = UniNorms, group = it, color = problem)) +
  theme_bw() +
  geom_line(aes(alpha = problem)) +
  scale_color_manual(values = c("black", "orange"), name = NULL) +
  scale_alpha_manual(values = c(0.1, 1)) +
  guides(alpha = "none", color = "none") +
  facet_wrap(vars(FPC), nrow = 4) +
  ggtitle("Univariate Norms of FPCs for Problematic Iterations")

