#' Family for Flexible Multivariate Joint Model
#'
#' This function specifies the different predictors and link functions as well
#' as the corresponding transform/updating/sampling functions as well as the
#' predict function.
#' @param ... All arguments are actually hard coded as needed by the
#'   implementation.
#'@export
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
    "optimizer" = MJM_opt,
    "sampler" = MJM_mcmc,
    "predict" = MJM_predict
  )
  
  class(rval) <- "family.bamlss"
  rval
}
