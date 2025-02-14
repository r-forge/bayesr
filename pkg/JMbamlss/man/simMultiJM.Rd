% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simMultiJM.R
\name{simMultiJM}
\alias{simMultiJM}
\title{New Simulation Function For Multivariate JMs Based On FPCs}
\usage{
simMultiJM(
  nsub = 300,
  times = seq(0, 120, 1),
  probmiss = 0.75,
  max_obs = length(times),
  maxfac = 1.5,
  nmark = 2,
  long_assoc = c("FPC", "splines", "param"),
  M = 6,
  FPC_bases = NULL,
  FPC_evals = NULL,
  mfpc_args = list(type = "split", eFunType = "Poly", ignoreDeg = NULL, eValType =
    "linear", eValScale = 1),
  re_cov_mat = NULL,
  ncovar = 2,
  lambda = function(t, x) {
     1.4 * log((t + 10)/1000)
 },
  gamma = function(x) {
     -1.5 + 0.3 * x[, 1]
 },
  alpha = rep(list(function(t, x) {
     0.3 + 0 * t
 }), nmark),
  mu = rep(list(function(t, x) {
     1.25 + 0.6 * sin(x[, 2]) + (-0.01) * t
 }),
    nmark),
  sigma = function(t, x) {
     0.3 + 0 * t + I(x$marker == "m2") * 0.2
 },
  tmax = NULL,
  seed = NULL,
  full = FALSE,
  file = NULL
)
}
\arguments{
\item{nsub}{Number of subjects.}

\item{times}{Vector of time points.}

\item{probmiss}{Probability of missingness.}

\item{max_obs}{Maximal number of observations per individual and marker. 
Defaults to no upper limit.}

\item{maxfac}{Factor changing the uniform censoring interval.}

\item{nmark}{Number of markers.}

\item{long_assoc}{Longitudinal association between the markers (Defaults to 
"FPC"). If "splines" or "param", then specify the normal covariance matrix
with argument 're_cov_mat' and include the random effects in argument mu.
 If "FPC", then principal components are used to model the association 
 structure.}

\item{M}{Number of principal components.}

\item{FPC_bases}{FunData object. If supplied, use the contained FPC as basis
for the association structure.}

\item{FPC_evals}{Vector of eigenvalues. If supplied, use the provided
eigenvalues for the association structure.}

\item{mfpc_args}{List containing the named arguments "type", "eFunType",
"ignoreDeg", "eValType" of function simMultiFunData and "eValScale" for
scaling the eigenvalues.}

\item{re_cov_mat}{If supplied, a covariance matrix to use for drawing the 
random effects needed for the association structure.}

\item{ncovar}{Number of covariates.}

\item{lambda}{Additive predictor of time-varying survival covariates.}

\item{gamma}{Additive predictor of time-constant survival covariates.}

\item{alpha}{List of length nmark containing the additive predictors of the
association.}

\item{mu}{List of length nmark containing the additive predictors of the 
longitudinal part.}

\item{sigma}{Additive predictor of the variance.}

\item{tmax}{Maximal time point of observations.}

\item{seed}{Seed for reproducibility.}

\item{full}{Create a wide-format data.frame and a short one containing only
survival info.}
}
\description{
Adapt the structure given by simJM function in bamlss.
}
