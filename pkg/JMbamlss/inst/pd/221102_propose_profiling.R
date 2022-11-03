
# PBC Data Example -------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/inst/")
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
source("R/preprocessing.R")
source("R/simMultiJM.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/MJM_opt.R")
source("R/opt_updating_cpp.R")
source("R/MJM_mcmc.R")
source("R/mcmc_proposing_cpp.R")
source("R/MJM_predict.R")
source("R/survint.R")
source("R/fpca.R")
source("R/mfpca_sim.R")
source("R/compile.R")
compile_alex(location)
sourceCpp("MatrixProd.cpp")

# .Rdata come from PBC analysis based on 221018_pbc_faster_update and are the
# objects for the function call of a propose_mjm in the first sampler after
# an optimizer of 15 iterations, with x being the first FPC term.
load(paste0(results_wd, "profile_mcmc.Rdata"))

library(lineprof)
l <- lineprof(propose_mjm(predictor, x, y, eta, eta_timegrid, eta_T, eta_T_mu,
                          eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_long,
                          eta_timegrid_lambda, survtime, logLik_old, nsubj, 
                          gq_weights, status, nmarker, nu, verbose_sampler,
                          prop))
shine(l)
# Apparently eigenMapMatMult is still the bottleneck of the computation, it is 
# found in the calculation of the Hesse matrix for t(X) %*% R^{-1} %*% X



library(microbenchmark)
delta <- rep(status, nmarker)
pred_l <- if(any(class(x) == "unc_pcre.random.effect")) "fpc_re" else "long"
int_i <- survint_C(pred = pred_l, pre_fac = exp(eta$gamma),
                   omega = exp(eta_timegrid),
                   int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                   weights = gq_weights, survtime = survtime)

set.seed(1702)
m <- microbenchmark(
  "score" = drop(crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
    crossprod(delta * x$XT, eta$alpha)) - int_i$score_int,
  "hess" = eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
    if (pred_l == "long" ) {
      matrix(int_i$hess_int, ncol = length(b_old))
    } else {
      diag(int_i$hess_int)
    }
)
boxplot(m, log = FALSE, ylab = "Milliseconds")

m_h <- microbenchmark(
  "Designmatrix" = eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X),
  "Survival Integral" = if (pred_l == "long" ) {
    matrix(int_i$hess_int, ncol = length(b_old))
  } else {
    diag(int_i$hess_int)
  }
)
boxplot(m_h, log = FALSE, ylab = "Microseconds")


mat <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X)
str(mat)
x$X[1:20, 1:8]
mat[1:20, 1:8]


library(Matrix)
m <- Matrix(x$X, sparse = TRUE)

# Sparse implementation very slow
m_s <- microbenchmark(
  "C++" = eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X),
  "Sparse" = crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
)

# Use psi_mat_crossprod function from:
#source("inst/pd/221103_develop_psi_mat_crossprod.R")
# Own R implementation very fast
m_o <- microbenchmark(
  "C++" = eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X),
  "R Loop" = psi_mat_crossprod(X = X, ni_obs = ni_obs, diags = diags)
)
m_o
