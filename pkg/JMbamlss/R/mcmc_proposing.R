
# Proposals for lambda ----------------------------------------------------

propose_mjm_lambda <- function(x, y, eta, eta_timegrid, eta_timegrid_lambda, 
                               eta_T, survtime, nu, ...){
  
  browser()
  b_old <- bamlss::get.state(x, "b")
  
  # Old logLikelihood and prior.
  gq_weights <- attr(y, "gq_weights")
  nsubj <- length(survtime)
  status <- attr(y, "status")

  sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
    (diag(nsubj)%x%t(gq_weights))%*%
    exp(eta_timegrid)
  logLik_old <- drop(status %*% eta_T - sum_Lambda) +
    sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
              log = TRUE))
  p_old <- x$prior(x$state$parameters)
  
  # Compute gradient and hessian integrals.
  int_i <- survint_gq(pred = "lambda", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_vec = x$Xgrid, weights = attr(y, "gq_weights"),
                      survtime = survtime)
  x_score <- drop(attr(y, "status") %*% x$XT) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  # Compute the inverse of the hessian.
  Sigma <- bamlss:::matrix_inv(x_H, index = NULL)
  
  # Get new location parameter for proposal
  mu_prop <- drop(b_old + nu * Sigma %*% x_score )
  
  # Sample new parameters.
  b_prop <- drop(rmvnorm(n = 1, mean = mu_prop, sigma = Sigma, method="chol"))
  names(b_prop) <- names(b_old)
  x$state$parameters <- bamlss::set.par(x$state$parameters, b_prop, "b")
  
  # Compute log priors.
  p_prop <- x$prior(x$state$parameters)
  qbeta_prop <- dmvnorm(b_prop, mean = mu_prop, sigma = Sigma, log = TRUE)
  
  # Update additive predictors
  fitted_timegrid_prop <- drop(x$Xgrid %*% b)
  eta_timegrid_lambda_prop <- eta_timegrid_lambda - x$state$fitted_timegrid + 
    fitted_timegrid_prop
  x$state$fitted_timegrid <- fitted_timegrid_prop
  eta_timegrid_prop <- eta_timegrid_lambda + eta_timegrid_long
  
  fit_prop <- drop(x$X %*% g)
  eta_prop <- eta
  eta_prop$lambda <- eta_prop$lambda - fitted(x$state) + fit_prop
  x$state$fitted.values <- fit_prop
  eta_T_prop <- eta_T - eta$lambda + eta_prop$lambda
  
  
  # New logLik
  sum_Lambda_prop <- (survtime/2 * exp(eta$gamma)) %*%
    (diag(nsubj)%x%t(gq_weights))%*%
    exp(eta_timegrid)
  logLik_prop <- drop(status %*% eta_T_prop - sum_Lambda_prop) +
    sum(dnorm(y[[1]][, "obs"], mean = eta_prop$mu, sd = exp(eta_prop$sigma),
              log = TRUE))
  
  
  # Compute gradient and hessian integrals for proposed beta
  int_i_prop <- survint_gq(pred = "lambda", pre_fac = exp(eta_prop$gamma),
                           omega = exp(eta_timegrid_prop),
                           int_vec = x$Xgrid, weights = attr(y, "gq_weights"),
                           survtime = survtime)
  x_score_prop <- drop(attr(y, "status") %*% x$XT) - 
    colSums(int_i_prop$score_int)
  x_H_prop <- matrix(colSums(int_i_prop$hess_int), ncol = b_prop)
  
  x_score_prop <- x_score_prop + x$grad(score = NULL, x$state$parameters,
                                        full = FALSE)
  x_H_prop <- x_H_prop + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  # Compute the inverse of the hessian.
  Sigma_prop <- bamlss:::matrix_inv(x_H, index = NULL)
  
  # Get new location parameter for proposal
  mu_prop <- drop(b_prop + nu * Sigma_prop %*% x_score_prop)
  
  
  # Prior prob
  qbeta <- dmvnorm(b_old, mean = mu_prop, sigma = Sigma_prop, log = TRUE)
  
  
  
  ## Save edf.
  x$state$edf <- sum_diag(int$hess %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  # cat(paste("\n",pibeta, pibetaprop, qbeta, qbetaprop, p1, p2, x$state$alpha, sep=";"))
  return(x$state)
}
