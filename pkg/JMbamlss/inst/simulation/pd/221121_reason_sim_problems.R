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
i <- 125
setting <- "scen_I_130922"
setwd("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,share=volkmana.hub/JMbamlss/")
# Load the data
load(paste0("simulation/", setting, "/data/d", i, ".Rdata"))
# Estimate the model using estimated FPCs
few_obs <- apply(table(d_rirs$data$id, d_rirs$data$marker), 1,
function (x) any(x < 2))
mfpca_est <- JMbamlss:::preproc_MFPCA(d_rirs$data %>%
filter(id %in% paste(which(!few_obs))) %>%
droplevels(),
uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
npc = 2, nbasis = 4)
vals <- which(mfpca_est$values > 0)
nfpc <- min(which(
cumsum(mfpca_est$values[vals])/sum(mfpca_est$values[vals]) > 0.95))
mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
list(functions = extractObs(mfpca$functions, i),
values = mfpca$values[i])
})
# Prepare objects for model fit
d_rirs_est <- JMbamlss:::attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)
f_est <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime, k = 20, bs = "ps"),
  gamma ~ 1 + x3,
  as.formula(paste0(
    "mu ~ -1 + marker + obstime:marker + x3:marker + obstime:x3:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_est_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)
# Model fit
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est,
timevar = "obstime", maxit = 1500, n.iter = 5500,
burnin = 500, thin = 5, verbose = TRUE)

b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est,
                timevar = "obstime", maxit = 466, sampler = FALSE,
                verbose = TRUE)
b_new <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est,
                timevar = "obstime", maxit = 500, sampler = FALSE, 
                start = parameters(b_est),
                verbose = TRUE)

# Censor survival time
d_rirs_est2 <- d_rirs_est %>%
  mutate(survtime = ifelse(survtime < 0.99, survtime, 0.99))
b_est <- bamlss(f_est, family = JMbamlss:::mjm_bamlss, data = d_rirs_est2,
                timevar = "obstime", maxit = 1500,  n.iter = 5500,
                burnin = 500, thin = 5)
