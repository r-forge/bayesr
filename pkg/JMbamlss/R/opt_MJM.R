

# Optimizer for MJM -------------------------------------------------------

opt_MJM <- function(x, y, start = NULL, eps = 0.0001, maxit = 100, nu = 0.1, 
                    ...) {
  
  if(!is.null(start))
    x <- bamlss:::set.starting.values(x, start)
  
  nmarker <- attr(y, "nmarker")
  marker <- attr(y, "marker")
  

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
      if(inherits(x$mu$smooth.construct[[j]], "pcre2.random.effect")){
        eta_grid_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$Xgrid))
      } else {
        eta_grid_sj <- drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
        eta_T_sj <- drop(x$mu$smooth.construct[[j]]$XT %*% b)
      }
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      x$mu$smooth.construct[[j]]$state$fitted_T <- eta_T_sj
      eta_timegrid_mu <- eta_timegrid_mu + eta_grid_sj
      eta_T_mu <- eta_T_mu + eta_T_sj
    }
  }
  
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      eta_sj <- drop(x$lambda$smooth.construct[[j]]$Xgrid %*% b)
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <- eta_sj
      eta_timegrid_lambda <- eta_timegrid_lambda + eta_sj
    }
  }
  
  if(!is.null(x$sigma$smooth.construct$model.matrix)) {
    x$sigma$smooth.construct$model.matrix$state$parameters[seq_len(nmarker)] <- 
      tapply(y[[1]][, "obs"], marker, function(x) log(sd(x, na.rm = TRUE)))
    x$sigma$smooth.construct$model.matrix$state$fitted.values <-
      x$sigma$smooth.construct$model.matrix$X %*% 
      x$sigma$smooth.construct$model.matrix$state$parameters
  }
  
  eta <- bamlss:::get.eta(x, expand = FALSE)
  eta_timegrid_long <- drop(
    t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
      (eta_timegrid_alpha*eta_timegrid_mu))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  
  # Für logLik
  nsubj <- attr(y, "nsubj")
  gq_weights <- attr(y, "gq_weights")
  take_last <- attr(y, "take_last")
  take_last_l <- attr(y, "take_last_l")
  status <- attr(y, "status")
  survtime <- y[[1]][, "time"][take_last]
  
  eps0 <- eps + 1
  eta0_surv <- do.call("cbind", eta[c("lambda", "gamma")])
  eta0_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
  eta0 <- cbind(eta0_surv, eta0_alpha)
  eta0_long <- do.call("cbind", eta[c("mu", "sigma")])
  iter <- 0
  
  # NUR ZUR NACHVERFOLGUNG DER GESCHÄTZTEN PARAMETER
  it_param <- list()
  
  # Updating the predictors
  while((eps0 > eps) & (iter < maxit)) {
    
    ## (1) update lambda.
    for(j in names(x$lambda$smooth.construct)) {
      state <- update_mjm_lambda(x$lambda$smooth.construct[[j]], y = y, nu = nu,
                                 eta = eta, eta_timegrid = eta_timegrid,
                                 survtime = survtime, ...)
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
                                  survtime = survtime, ...)
        eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[j]]$state) +
          fitted(state)
        x$gamma$smooth.construct[[j]]$state <- state
      }
    }
    
    ## (3) update alpha.
    if(length(x$alpha$smooth.construct)) {
      for(j in seq_along(x$alpha$smooth.construct)) {
        state <- update_mjm_alpha(x$alpha$smooth.construct[[j]], y = y, nu = nu,
                                  eta = eta, eta_timegrid = eta_timegrid,
                                  eta_timegrid_mu = eta_timegrid_mu,
                                  eta_T_mu = eta_T_mu, survtime = survtime, ...)
        eta$alpha <- eta$alpha -
          drop(fitted(x$alpha$smooth.construct[[j]]$state)) +
          fitted(state)
        eta_timegrid_alpha <- eta_timegrid_alpha -
          x$alpha$smooth.construct[[j]]$state$fitted_timegrid +
          state$fitted_timegrid
        eta_timegrid_long <- drop(
          t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
            (eta_timegrid_alpha * eta_timegrid_mu))
        eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
        x$alpha$smooth.construct[[j]]$state <- state
      }
    }
    
    ## (4) update mu.
    if(length(x$mu$smooth.construct)) {
      for(j in seq_along(x$mu$smooth.construct)) {
        state <- update_mjm_mu(x$mu$smooth.construct[[j]], y = y, nu = nu,
                               eta = eta, eta_timegrid = eta_timegrid,
                               eta_timegrid_alpha = eta_timegrid_alpha,
                               survtime = survtime, ...)
        eta$mu <- eta$mu -
          drop(fitted(x$mu$smooth.construct[[j]]$state)) +
          fitted(state)
        eta_timegrid_mu <- eta_timegrid_mu -
          x$mu$smooth.construct[[j]]$state$fitted_timegrid +
          state$fitted_timegrid
        eta_timegrid_long <- drop(
          t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
            (eta_timegrid_alpha * eta_timegrid_mu))
        eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
        eta_T_mu <- eta_T_mu -
          x$mu$smooth.construct[[j]]$state$fitted_T +
          state$fitted_T
        x$mu$smooth.construct[[j]]$state <- state
      }
    }
    
    # Likelihood calculation
    # Was passiert, wenn es keine longitudinale Beobachtung gibt für den Event-
    # Zeitpunkt? Hier bräuchte man eigentlich alpha und mu als nsubj*nmarker
    # Vektor.
    eta_T_long <- drop(
      t(rep(1, nmarker)) %x% diag(nsubj) %*% (eta$alpha*eta_T_mu))
    eta_T <- eta$lambda + eta$gamma + eta_T_long
    
    sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
      (diag(nsubj)%x%t(gq_weights))%*%
      exp(eta_timegrid)
    logLik <- drop(status %*% eta_T - sum_Lambda) +
      sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
                log = TRUE))
    # Eigentlich sollte hier doch auch über die Log-Posterior das Max gebildet
    # werden? Die Score und Hesse-Funktionen beziehen nämlich schon die Prioris
    # mit ein
    
    # Stopping criterion
    iter <- iter + 1
    eta1_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eta1_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
    eta1 <- cbind(eta1_surv, eta0_alpha)
    eps0_surv <- mean(abs((eta1 - eta0) / eta1), na.rm = TRUE)
    eta1_long <- do.call("cbind", eta[c("mu", "sigma")])
    eps0_long <- mean(abs((eta1_long - eta0_long) / eta1_long), na.rm = TRUE)
    eps0 <- mean(c(eps0_surv, eps0_long))
    
    #if (iter %% 5 == 0) {
    cat("It ", iter,", LogLik ", logLik, "\n")
    #}
    eta0 <- eta1
    eta0_long <- eta1_long
    
    # NUR ZUR NACHVERFOLGUNG DER GESCHÄTZTEN PARAMETER
    it_param[[iter]] <- bamlss:::get.all.par(x)
  }
  
  # Log-Posterior ausrechnen und ausgeben
  # get_logPost Funktion aus JM als Vorlage
  ## return(list("parameters" = par, "fitted.values" = eta))
  
  # max_y <- max(y[[1]][, 1])
  # pred_data <- data.frame(seq(0, max_y, length.out = 100))
  # colnames(pred_data) <- x$lambda$smooth.construct[[1]]$term
  # k <- x$lambda$smooth.construct[[1]]$bs.dim
  # b_it <- x$lambda$smooth.construct[[1]]$state$parameters[-k]
  # pred_mat <- PredictMat(x$lambda$smooth.construct[[1]], pred_data)
  # plot(pred_data[, 1], pred_mat%*%b_it)
  
  assign("it_param", it_param, envir = .GlobalEnv)
  logPost <- as.numeric(logLik + bamlss:::get.log.prior(x))
  return(list("fitted.values" = eta, "parameters" = bamlss:::get.all.par(x),
              "logLik" = logLik, "logPost" = logPost,
              "hessian" = bamlss:::get.hessian(x),
              "converged" = iter < maxit))
  
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