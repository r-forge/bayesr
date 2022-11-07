
# Proposals for all predictors --------------------------------------------

propose_mjm_old <- function(predictor, x, y, eta, eta_timegrid, eta_T, eta_T_mu,
                        eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_long,
                        eta_timegrid_lambda, survtime, logLik_old, nsubj, 
                        gq_weights, status, nmarker, nu, verbose_sampler,
                        prop#, prop_dens
                        ) {
  
  #if (predictor == "alpha") {
    # # Old logLik -- Why does it change so drastically?
    # sum_Lambda_old <- (survtime/2 * exp(eta$gamma)) %*%
    #   (diag(nsubj)%x%t(gq_weights))%*%
    #   exp(eta_timegrid)
    # logLik_old_comp <- drop(status %*% eta_T - sum_Lambda_old) +
    #   sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
    #             log = TRUE))
    # eta_old <- eta
    # eta_timegrid_old <- eta_timegrid
    # eta_T_old <- eta_T
    # 
    # # Apparently eta_timegrid has changed drastically, so sum_Lambda is now 
    # # crazy high
    # state_fitted_timegrid_old <- x$state$fitted_timegrid
    # 
  #}
  n_w <- length(gq_weights)
  
  # Get old parameters
  b_old <- bamlss::get.state(x, "b")
  p_old <- x$prior(x$state$parameters)
  
  if (predictor %in% c("alpha", "mu")) {
    delta <- rep(status, nmarker)
  }
  
  if (predictor == "mu") {
    pred_l <- if(any(class(x) == "unc_pcre.random.effect")) "fpc_re" else "long"
  }
  
  # Calculate Newton step based on old parameter
  #if (is.null(prop_dens)) {
    switch(predictor,
           "lambda" = {
             int_i <- survint_C(pred = "lambda", pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_vec = x$Xgrid, weights = gq_weights,
                                survtime = survtime)
             x_score <- drop(status %*% x$XT) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "gamma" = {
             int_i <- survint_C(pred = "gamma", pre_fac = exp(eta$gamma), 
                                pre_vec = x$X, omega = exp(eta_timegrid),
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(status %*% x$X) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "alpha" = {
             int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(t(delta * x$XT) %*% eta_T_mu) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "mu" = {
             int_i <- survint_C(pred = pred_l, pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(
               crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
                 crossprod(delta * x$XT, eta$alpha)) - int_i$score_int
             x_H <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
               if (pred_l == "long" ) {
                 matrix(int_i$hess_int, ncol = length(b_old))
               } else {
                 diag(int_i$hess_int)
               }
           },
           "sigma" = {
             x_score <- crossprod(x$X, -1 + (y[[1]][, "obs"] - eta$mu)^2 / 
                                    exp(eta$sigma)^2)
             x_H <- 2 * crossprod(x$X * 
                                    drop((y[[1]][, "obs"] - eta$mu) / 
                                           exp(eta$sigma)^2),
                                  x$X * drop(y[[1]][, "obs"] - eta$mu))
           })
    
    # Include priori in derivatives
    x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
    x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
    
    # Compute the inverse of the hessian.
    Sigma_prop <- bamlss:::matrix_inv(x_H, index = NULL)
    
    # Get new location parameter for proposal
    mu_prop <- drop(b_old + nu * Sigma_prop %*% x_score )
    
  # } else {
  #   
  #   mu_prop <- prop_dens[[1]]$mu
  #   Sigma_prop <- prop_dens[[1]]$Sigma
  #   
  # }
  
  
  
  # Sample new parameters.
  # if(is.null(prop)) {
    b_prop <- drop(mvtnorm::rmvnorm(n = 1, mean = mu_prop, sigma = Sigma_prop,
                                    method="chol"))
  # } else {
  #   b_prop <- prop[seq_along(mu_prop)]
  # }

  names(b_prop) <- names(b_old)
  x$state$parameters <- bamlss::set.par(x$state$parameters, b_prop, "b")
  p_prop <- x$prior(x$state$parameters)
  
  # Compute proposal density of proposed parameter given old
  q_prop_giv_old <- mvtnorm::dmvnorm(b_prop, mean = mu_prop, sigma = Sigma_prop,
                                     log = TRUE)
  
  
  # Update additive predictors and x$state
  switch(predictor,
    "lambda" = {
      
      # timegrids
      fitted_timegrid_prop <- drop(x$Xgrid %*% b_prop)
      eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + 
        fitted_timegrid_prop
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      
      # fitted values
      fit_prop <- drop(x$X %*% b_prop)
      eta_T <- eta_T - eta$lambda
      eta$lambda <- eta$lambda - fitted(x$state) + fit_prop
      eta_T <- eta_T + eta$lambda
      
      # state
      x$state$fitted_timegrid <- fitted_timegrid_prop
      x$state$fitted.values <- fit_prop
      
    },
    "gamma" = {
      
      # fitted values and state
      fit_prop <- drop(x$X %*% b_prop)
      eta_T <- eta_T - eta$gamma
      eta$gamma <- eta$gamma - fitted(x$state) + fit_prop
      eta_T <- eta_T + eta$gamma
      x$state$fitted.values <- fit_prop
      
    },
    "alpha" = {
      
      # timegrids
      fitted_timegrid_prop <- drop(x$Xgrid %*% b_prop)
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + 
        fitted_timegrid_prop
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu, 
                                          nrow = nsubj*n_w, ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
      
      # fitted values
      fit_prop <- drop(x$X %*% b_prop)
      eta$alpha <- eta$alpha - fitted(x$state) + fit_prop
      eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj, 
                                   ncol = nmarker))
      eta_T <- eta$lambda + eta$gamma + eta_T_long
      
      # state
      x$state$fitted_timegrid <- fitted_timegrid_prop
      x$state$fitted.values <- fit_prop
      
    },
    "mu" = {
     
      # timegrids
      fitted_timegrid_prop <- drop(x$Xgrid %*% b_prop)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + 
       fitted_timegrid_prop
      eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu,
                                          nrow = nsubj*n_w, ncol = nmarker))
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
     
      # fitted longitudinal values
      fit_prop <- drop(x$X %*% b_prop)
      eta$mu <- eta$mu - fitted(x$state) + fit_prop
      
      # fitted survival values
      fit_T_prop <- drop(x$XT %*% b_prop)
      eta_T_mu <- eta_T_mu - x$state$fitted_T + fit_T_prop
      eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj, 
                                   ncol = nmarker))
      eta_T <- eta$lambda + eta$gamma + eta_T_long
     
     # state
     x$state$fitted_timegrid <- fitted_timegrid_prop
     x$state$fitted.values <- fit_prop
     x$state$fitted_T <- fit_T_prop
     
    },
    "sigma" = {
      
      # fitted values and state
      fit_prop <- drop(x$X %*% b_prop)
      eta$sigma <- eta$sigma - fitted(x$state) + fit_prop
      x$state$fitted.values <- fit_prop
      
    })
  
  # New logLik
  sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
    (diag(nsubj)%x%t(gq_weights))%*%
    exp(eta_timegrid)
  logLik <- drop(status %*% eta_T - sum_Lambda) +
    sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
              log = TRUE))
  
  # Calculate Newton step based on proposed parameter and updated etas
  # if (is.null(prop_dens)) {
  #   
    switch(predictor,
           "lambda" = {
             int_i <- survint_C(pred = "lambda", pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_vec = x$Xgrid, weights = gq_weights,
                                survtime = survtime)
             x_score <- drop(status %*% x$XT) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "gamma" = {
             int_i <- survint_C(pred = "gamma", pre_fac = exp(eta$gamma), 
                                pre_vec = x$X, omega = exp(eta_timegrid),
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(status %*% x$X) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "alpha" = {
             int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(t(delta * x$XT) %*% eta_T_mu) - int_i$score_int
             x_H <- matrix(int_i$hess_int, ncol = length(b_old))
           },
           "mu" = {
             int_i <- survint_C(pred = pred_l, pre_fac = exp(eta$gamma),
                                omega = exp(eta_timegrid),
                                int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                                weights = gq_weights, survtime = survtime)
             x_score <- drop(
               crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
                 crossprod(delta * x$XT, eta$alpha)) - int_i$score_int
             x_H <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
               if (pred_l == "long" ) {
                 matrix(int_i$hess_int, ncol = length(b_old))
               } else {
                 diag(int_i$hess_int)
               }
           },
           "sigma" = {
             x_score <- crossprod(x$X, -1 + (y[[1]][, "obs"] - eta$mu)^2 / 
                                    exp(eta$sigma)^2)
             x_H <- 2 * crossprod(x$X * 
                                    drop((y[[1]][, "obs"] - eta$mu) / 
                                           exp(eta$sigma)^2),
                                  x$X * drop(y[[1]][, "obs"] - eta$mu))
             x_H0 <- x_H
           })
    x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
    x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
    Sigma <- bamlss:::matrix_inv(x_H, index = NULL)
    
    # Get new location parameter for proposal
    mu <- drop(b_prop + nu * Sigma %*% x_score)
    
  # } else {
  #   
  #   mu <- prop_dens[[2]]$mu
  #   Sigma <- prop_dens[[2]]$Sigma
  #   
  #   if (predictor == "sigma") {
  #     x_H0 <- prop_dens[[2]]$edf
  #   } else {
  #     int_i <- prop_dens[[2]]$edf
  #   }
  #   
  # }
  
  
  # Proposal density of old parameter given proposed
  q_old_giv_prop <- mvtnorm::dmvnorm(b_old, mean = mu, sigma = Sigma, 
                                     log = TRUE)
  
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((logLik + q_old_giv_prop + p_prop) - 
                          (logLik_old + q_prop_giv_old + p_old))
  
  
  # WIE KOMMTS?
  ## Save edf.
  x$state$edf <- if (predictor == "sigma") {
    bamlss:::sum_diag(x_H0 %*% Sigma)
  } else if (predictor == "mu") {
    if (pred_l == "fpc_re") {
      bamlss:::sum_diag(diag(int_i$hess_int) %*% Sigma)
    } else {
      bamlss:::sum_diag(matrix(int_i$hess_int, ncol = length(b_prop)) %*%
                          Sigma)
    }
  } else {
    bamlss:::sum_diag(matrix(int_i$hess_int, ncol = length(b_prop)) %*%
                        Sigma)
  }
  
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if (!is.null(prop)) {
      taus <- grep("tau2", names(prop))
      x$state$parameters <- bamlss::set.par(x$state$parameters, prop[taus], 
                                            "tau2")
    } else {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(b_prop, x$S[[1]]) %*% b_prop + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- bamlss::set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- bamlss:::uni.slice(
          x$state$parameters, x, NULL, NULL, NULL, id = predictor, j, 
          logPost = bamlss:::uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    } }
  }
  
  if(verbose_sampler) {
    cat(predictor, "LLO:", logLik_old, "LLN:", logLik, 
        "PropO:", q_old_giv_prop, "PropN:", q_prop_giv_old, 
        "PriO:", p_old, "PriN:", p_prop, "Alph:", exp(x$state$alpha), "\n")
  }
  
  return(list(xstate = x$state,
              etas = switch(predictor,
                            "lambda" = list(eta = eta, eta_T = eta_T, 
                                            eta_timegrid = eta_timegrid, 
                                            eta_timegrid_lambda = 
                                              eta_timegrid_lambda),
                            "gamma" = list(eta = eta, eta_T = eta_T),
                            "alpha" = list(eta = eta, eta_T = eta_T, 
                                           eta_timegrid = eta_timegrid,
                                           eta_timegrid_long = 
                                             eta_timegrid_long,
                                           eta_timegrid_alpha = 
                                             eta_timegrid_alpha),
                            "mu" = list(eta = eta, eta_T = eta_T, 
                                        eta_T_mu = eta_T_mu,
                                        eta_timegrid = eta_timegrid,
                                        eta_timegrid_long = 
                                          eta_timegrid_long,
                                        eta_timegrid_mu = eta_timegrid_mu),
                            "sigma" = list(eta = eta)),
              logLik = logLik
              # ,prop_dens = list(list("mu" = mu_prop, "Sigma" = Sigma_prop),
              #                   list("mu" = mu, "Sigma" = Sigma,
              #                        "edf" = if (predictor == "sigma") {
              #                          x_H0} else int_i))
              ))
  
}



# Proposals for lambda ----------------------------------------------------

propose_mjm_lambda <- function(x, y, eta, eta_timegrid, eta_timegrid_lambda, 
                               eta_T, survtime, nu, eta_timegrid_long, ...){
  
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
  x_score <- drop(attr(y, "status") %*% x$XT) - int_i$score_int
  x_H <- matrix(int_i$hess_int, ncol = length(b_old))
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  # Compute the inverse of the hessian.
  Sigma <- bamlss:::matrix_inv(x_H, index = NULL)
  
  # Get new location parameter for proposal
  mu_prop <- drop(b_old + nu * Sigma %*% x_score )
  
  # Sample new parameters.
  b_prop <- drop(mvtnorm::rmvnorm(n = 1, mean = mu_prop, sigma = Sigma,
                                  method="chol"))
  names(b_prop) <- names(b_old)
  x$state$parameters <- bamlss::set.par(x$state$parameters, b_prop, "b")
  
  # Compute log priors.
  p_prop <- x$prior(x$state$parameters)
  q_prop_giv_old <- mvtnorm::dmvnorm(b_prop, mean = mu_prop, sigma = Sigma, 
                                 log = TRUE)
  
  # Update additive predictors
  fitted_timegrid_prop <- drop(x$Xgrid %*% b_prop)
  eta_timegrid_lambda_prop <- eta_timegrid_lambda - x$state$fitted_timegrid + 
    fitted_timegrid_prop
  x$state$fitted_timegrid <- fitted_timegrid_prop
  eta_timegrid_prop <- eta_timegrid_lambda + eta_timegrid_long
  
  fit_prop <- drop(x$X %*% b_prop)
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
    int_i_prop$score_int
  x_H_prop <- matrix(int_i_prop$hess_int, ncol = length(b_prop))
  
  x_score_prop <- x_score_prop + x$grad(score = NULL, x$state$parameters,
                                        full = FALSE)
  x_H_prop <- x_H_prop + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  # Compute the inverse of the hessian.
  Sigma_prop <- bamlss:::matrix_inv(x_H_prop, index = NULL)
  
  # Get new location parameter for proposal
  mu_prop <- drop(b_prop + nu * Sigma_prop %*% x_score_prop)
  
  
  # Prior prob
  q_old_giv_prop <- mvtnorm::dmvnorm(b_old, mean = mu_prop, sigma = Sigma_prop, 
                                log = TRUE)
  
  ## Save edf.
  x$state$edf <- bamlss:::sum_diag(
    matrix(int_i_prop$hess_int, ncol = length(b_prop)) %*% Sigma_prop)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- bamlss::get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- bamlss::set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- bamlss:::uni.slice(
          x$state$parameters, x, NULL, NULL, NULL, id = "mu", j, 
          logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((logLik_prop + q_old_giv_prop + p_prop) - 
                          (logLik_old + q_prop_giv_old + p_old))
  
  return(x$state)
}
