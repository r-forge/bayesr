

# Optimizer for MJM -------------------------------------------------------

MJM_opt <- function(x, y, start = NULL, eps = 0.0001, maxit = 100, nu = 0.1,
                    opt_long = TRUE, alpha.eps = 0.001, par_trace = FALSE,
                    verbose = FALSE, update_nu = FALSE, update_tau = FALSE,
                    ...) {
  
  if(!is.null(start))
    x <- bamlss:::set.starting.values(x, start)
  
  nsubj <- attr(y, "nsubj")
  gq_weights <- attr(y, "gq_weights")
  n_w <- length(gq_weights)
  take_last <- attr(y, "take_last")
  take_last_l <- attr(y, "take_last_l")
  status <- attr(y, "status")
  survtime <- y[[1]][, "time"][take_last]
  nmarker <- attr(y, "nmarker")
  marker <- attr(y, "marker")
  logLik <- NULL
  edf <- 0
  

  ## Set alpha/mu/sigma intercept starting value.
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$alpha$smooth.construct[[j]]$state$parameters[seq_len(nmarker)] <- 
          .Machine$double.eps
        x$alpha$smooth.construct[[j]]$state$fitted.values <- 
          x$alpha$smooth.construct[[j]]$X %*% 
          x$alpha$smooth.construct[[j]]$state$parameters
      } 
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      eta_grid_sj <- drop(x$alpha$smooth.construct[[j]]$Xgrid %*% b)
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      eta_timegrid_alpha <- eta_timegrid_alpha + eta_grid_sj
      edf <- edf + x$alpha$smooth.construct[[j]]$state$edf
    }
  } 
  eta_timegrid_mu <- 0
  eta_T_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$mu$smooth.construct[[j]]$state$parameters[seq_len(nmarker)] <- 
          tapply(y[[1]][, "obs"], marker, mean, na.rm = TRUE)
        x$mu$smooth.construct[[j]]$state$fitted.values <- 
          x$mu$smooth.construct[[j]]$X %*% 
          x$mu$smooth.construct[[j]]$state$parameters
      }
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      ### PASST DAS HIER NOCH MIT DEM PCRE2 ?
      if(inherits(x$mu$smooth.construct[[j]], "pcre2.random.effect")){
        eta_grid_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$Xgrid))
        eta_T_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$XT))
      } else {
        eta_grid_sj <- drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
        eta_T_sj <- drop(x$mu$smooth.construct[[j]]$XT %*% b)
      }
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      x$mu$smooth.construct[[j]]$state$fitted_T <- eta_T_sj
      eta_timegrid_mu <- eta_timegrid_mu + eta_grid_sj
      eta_T_mu <- eta_T_mu + eta_T_sj
      edf <- edf + x$mu$smooth.construct[[j]]$state$edf
    }
  }
  
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      eta_sj <- drop(x$lambda$smooth.construct[[j]]$Xgrid %*% b)
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <- eta_sj
      eta_timegrid_lambda <- eta_timegrid_lambda + eta_sj
      edf <- edf + x$lambda$smooth.construct[[j]]$state$edf
    }
  }
  
  if(!is.null(x$sigma$smooth.construct$model.matrix)) {
    if(is.null(start)) {
      x$sigma$smooth.construct$model.matrix$state$parameters[
        seq_len(nmarker)] <- tapply(y[[1]][, "obs"], marker,
                                    function(x) log(sd(x, na.rm = TRUE)))
    }
    x$sigma$smooth.construct$model.matrix$state$fitted.values <-
      x$sigma$smooth.construct$model.matrix$X %*% 
      x$sigma$smooth.construct$model.matrix$state$parameters
    edf <- edf + x$sigma$smooth.construct$model.matrix$state$edf
  }
  
  if (length(x$gamma$smooth.construct)) {
    for (j in names(x$gamma$smooth.construct)) {
      edf <- edf + x$gamma$smooth.construct[[j]]$state$edf
    } 
  }
  
  eta <- bamlss:::get.eta(x, expand = FALSE)
  eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu,
                                      nrow = nsubj*n_w, ncol = nmarker))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  
  
  ## Extract current value of the log-posterior.
  # Define here so data is available in the function environment
  # Used for updating nu / verbose 
  get_LogLik <- function(eta_timegrid, eta_T_mu, eta) {
    
    eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj,
                                 ncol = nmarker))
    eta_T <- eta$lambda + eta$gamma + eta_T_long
    sum_Lambda <- drop(
      crossprod(survtime/2 * exp(eta$gamma), 
                colSums(gq_weights*matrix(exp(eta_timegrid), 
                                          ncol = nsubj, 
                                          nrow = length(gq_weights)))))
    logLik <- drop(crossprod(status, eta_T)) - sum_Lambda +
      sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
                log = TRUE))
    
    return(logLik)

  }
  
  # For algorithm
  eps0 <- eps0_surv <- eps0_long <- eps + 1
  eta0_surv <- do.call("cbind", eta[c("lambda", "gamma")])
  eta0_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
  eta0 <- cbind(eta0_surv, eta0_alpha)
  eta0_long <- do.call("cbind", eta[c("mu", "sigma")])
  iter <- 0
  alpha_update <- if(is.null(start)) FALSE else TRUE
  
  # NUR ZUR NACHVERFOLGUNG DER GESCHÄTZTEN PARAMETER
  if (par_trace) {
    it_param <- list()
  }
  
  # Updating the predictors
  while((eps0 > eps) & (iter < maxit)) {
   
    ## (1) update lambda.
    for(j in names(x$lambda$smooth.construct)) {
      state <- update_mjm_lambda(x$lambda$smooth.construct[[j]], y = y, nu = nu,
                                 eta = eta, eta_timegrid = eta_timegrid,
                                 eta_T_mu = eta_T_mu,
                                 survtime = survtime, update_nu = update_nu,
                                 get_LogLik = get_LogLik,
                                 update_tau = update_tau, edf = edf, ...)
      eta_timegrid_lambda <- eta_timegrid_lambda -
        x$lambda$smooth.construct[[j]]$state$fitted_timegrid +
        state$fitted_timegrid
      eta_timegrid <- eta_timegrid_lambda +
        eta_timegrid_long
      eta$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[j]]$state) +
        fitted(state)
      x$lambda$smooth.construct[[j]]$state <- state
    }
    
    ## (2) update gamma.
    if(length(x$gamma$smooth.construct)) {
      for(j in seq_along(x$gamma$smooth.construct)) {
        state <- update_mjm_gamma(x$gamma$smooth.construct[[j]], y = y, nu = nu,
                                  eta = eta, eta_timegrid = eta_timegrid,
                                  eta_T_mu = eta_T_mu,
                                  survtime = survtime, update_nu = update_nu,
                                  get_LogLik = get_LogLik,
                                  update_tau = update_tau, edf = edf, ...)
        eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[j]]$state) +
          fitted(state)
        x$gamma$smooth.construct[[j]]$state <- state
      }
    }
    
    if (!alpha_update && max(c(eps0_surv, eps0_long)) < alpha.eps) {
      alpha_update <- TRUE
      if (verbose) {
        cat("It ", iter, "-- Start Alpha Update", "\n") 
      }
    }
    if (opt_long) {
      ## (3) update alpha.
      if(alpha_update) {
        if(length(x$alpha$smooth.construct)) {
          for(j in seq_along(x$alpha$smooth.construct)) {
            state <- update_mjm_alpha(x$alpha$smooth.construct[[j]], y = y, 
                                      nu = nu,
                                      eta = eta, eta_timegrid = eta_timegrid,
                                      eta_timegrid_lambda = eta_timegrid_lambda,
                                      eta_timegrid_mu = eta_timegrid_mu,
                                      eta_timegrid_alpha = eta_timegrid_alpha,
                                      eta_T_mu = eta_T_mu, survtime = survtime, 
                                      update_nu = update_nu,
                                      get_LogLik = get_LogLik,
                                      update_tau = update_tau, edf = edf, ...)
            eta$alpha <- eta$alpha -
              drop(fitted(x$alpha$smooth.construct[[j]]$state)) +
              fitted(state)
            eta_timegrid_alpha <- eta_timegrid_alpha -
              x$alpha$smooth.construct[[j]]$state$fitted_timegrid +
              state$fitted_timegrid
            eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha *
                                                  eta_timegrid_mu,
                                                nrow = nsubj*n_w,
                                                ncol = nmarker))
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
            x$alpha$smooth.construct[[j]]$state <- state
          }
        }
      }
      
      ## (4) update mu.
      if(length(x$mu$smooth.construct)) {
        for(j in seq_along(x$mu$smooth.construct)) {
          state <- update_mjm_mu(x$mu$smooth.construct[[j]], y = y, nu = nu,
                                 eta = eta, eta_timegrid = eta_timegrid,
                                 eta_timegrid_lambda = eta_timegrid_lambda,
                                 eta_timegrid_mu = eta_timegrid_mu,
                                 eta_timegrid_alpha = eta_timegrid_alpha,
                                 eta_T_mu = eta_T_mu, survtime = survtime,
                                 get_LogLik = get_LogLik, update_nu = update_nu,
                                 update_tau = update_tau, edf = edf, ...)
          eta$mu <- eta$mu -
            drop(fitted(x$mu$smooth.construct[[j]]$state)) + fitted(state)
          eta_timegrid_mu <- eta_timegrid_mu -
            x$mu$smooth.construct[[j]]$state$fitted_timegrid +
            state$fitted_timegrid
          eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha * 
                                                eta_timegrid_mu, 
                                              nrow = nsubj*n_w, ncol = nmarker))
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
          eta_T_mu <- eta_T_mu -
            x$mu$smooth.construct[[j]]$state$fitted_T +
            state$fitted_T
          x$mu$smooth.construct[[j]]$state <- state
        }
      }
      
      ## (5) update sigma.
      if(length(x$sigma$smooth.construct)) {
        for(j in seq_along(x$sigma$smooth.construct)) {
          state <- update_mjm_sigma(x$sigma$smooth.construct[[j]], y = y, 
                                    nu = nu, eta = eta, 
                                    eta_timegrid = eta_timegrid,
                                    eta_T_mu = eta_T_mu,
                                    get_LogLik = get_LogLik,
                                    survtime = survtime, update_nu = update_nu,
                                    update_tau = update_tau, edf = edf, ...)
          eta$sigma <- eta$sigma -
            drop(fitted(x$sigma$smooth.construct[[j]]$state)) +
            fitted(state)
          x$sigma$smooth.construct[[j]]$state <- state
        }
      }
    }

   
    if (verbose) {
      # Likelihood calculation
      # Was passiert, wenn es keine longitudinale Beobachtung gibt für den Event-
      # Zeitpunkt? Hier bräuchte man eigentlich alpha und mu als nsubj*nmarker
      # Vektor.
      
      logLik <- get_LogLik(eta_timegrid = eta_timegrid,
                           eta_T_mu = eta_T_mu,
                           eta = eta)
      
      cat("It ", iter,", LogLik ", logLik, ", Post", 
          as.numeric(logLik + bamlss:::get.log.prior(x)), "\n")
    }
    
    
    # Stopping criterion
    iter <- iter + 1
    eta1_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eta1_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
    eta1_long <- do.call("cbind", eta[c("mu", "sigma")])
    eps0_surv <- mean(abs((eta1_surv - eta0_surv) / eta1_surv), na.rm = TRUE)
    eps0_alpha <- if (alpha_update) {
      mean(abs((eta1_alpha - eta0_alpha) / eta1_alpha), na.rm = TRUE)
    } else {
      eps + 1
    }
    eps0_long <- mean(abs((eta1_long - eta0_long) / eta1_long), na.rm = TRUE)
    eps0 <- mean(c(eps0_surv, eps0_alpha, eps0_long))
    
    eta0_surv <- eta1_surv
    eta0_alpha <- eta1_alpha
    eta0_long <- eta1_long
    
    # NUR ZUR NACHVERFOLGUNG DER GESCHÄTZTEN PARAMETER
    if (par_trace) {
      it_param[[iter]] <- bamlss:::get.all.par(x)
    }
    
  }
  
  # Log-Posterior ausrechnen und ausgeben
  if (is.null(logLik)) {
    logLik <- get_LogLik(eta_timegrid = eta_timegrid,
                         eta_T_mu = eta_T_mu,
                         eta = eta)
  }
  logPost <- as.numeric(logLik + bamlss:::get.log.prior(x))
  return(list("fitted.values" = eta, "parameters" = bamlss:::get.all.par(x),
              "logLik" = logLik, "logPost" = logPost,
              "hessian" = bamlss:::get.hessian(x),
              "converged" = iter < maxit, 
              "par_trace" = if (par_trace) it_param else NULL))
  
  assign("b", x, envir = .GlobalEnv)
  assign("eta", list(eta = eta,
                     eta_timegrid = eta_timegrid,
                     eta_timegrid_long = eta_timegrid_long,
                     eta_timegrid_alpha = eta_timegrid_alpha, 
                     eta_timegrid_mu = eta_timegrid_mu,
                     eta_timegrid_lambda = eta_timegrid_lambda,
                     eta_T_mu = eta_T_mu),
         envir = .GlobalEnv)
  stop("Basst.")
  
}
