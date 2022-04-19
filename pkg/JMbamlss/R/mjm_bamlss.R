mjm_bamlss <- function(...)
{
  links = c(
    lambda = "log",
    gamma = "log",
    mu = "identity",
    sigma = "log",
    alpha = "identity"
  )
  
  rval <- list(
    "family" = "mjm",
    "names" = c("lambda", "gamma", "mu", "sigma", "alpha"),
    "links" = links,
    "transform" = MJM_transform,
    "optimizer" = opt_MJM,
    "sampler" = MJM_mcmc,
    "predict" = MJM_predict
  )
  
  class(rval) <- "family.bamlss"
  rval
}
