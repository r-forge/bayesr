
# MJM_mcmc ----------------------------------------------------------------

MJM_mcmc <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
                     n.iter = 1200, burnin = 200, thin = 1, step = 20, 
                     nu_sampler = 1, prop_pred = NULL, verbose_sampler = FALSE,
                     prop_list = NULL, ...)
{
########## REMOVE prop_pred
  # Set starting values for the sampling
  if(!is.null(start)) {
    if(is.matrix(start)) {
      if(any(i <- grepl("Mean", colnames(start))))
        start <- start[, i]
      else stop("the starting values should be a vector not a matrix!")
    }
    x <- bamlss:::set.starting.values(x, start)
  }
  
  # Names of parameters/predictors and other attributes
  nx <- if (is.null(prop_pred)) names(x) else prop_pred
  nmarker <- attr(y, "nmarker")
  take_last <- attr(y, "take_last")
  survtime <- y[[1]][, "time"][take_last]
  nsubj <- length(survtime)
  gq_weights <- attr(y, "gq_weights")
  nw <- length(gq_weights)
  status <- attr(y, "status")
  
  ## Number of observations.
  #nobs <- attr(y, "nobs")
  
  ## Number of subdivions used for the time grid.
  #sub <- attr(y, "subdivisions")
  
  ## The interval width from subdivisons.
  #width <- attr(y, "width")
  
  ## Subject specific indicator
  # take <- attr(y, "take")
  # nlong <- length(take)
  
  ## Extract the status for individual i.
  # status <- y[take, "status"]
  
  # ## nonlinear setup
  # nonlinear <- attr(y, "nonlinear")
  # tp <- attr(y, "tp") 
  # fac <- attr(y, "fac")
  
  # ## Make id for individual i.
  # id <- which(take)
  # id <- append(id[-1], nlong + 1) - id 
  # id <- rep(1:nobs, id)
  
  ## Compute additive predictors.
  eta <- bamlss:::get.eta(x, expand = FALSE)
  
  
  ## For the time dependent part, compute
  ## predictors based on the time grid.
  # mu
  eta_timegrid_mu <- 0
  eta_T_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- 
        drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
      x$mu$smooth.construct[[j]]$state$fitted_T <- 
        drop(x$mu$smooth.construct[[j]]$XT %*% b)
      eta_timegrid_mu <- eta_timegrid_mu + 
        x$mu$smooth.construct[[j]]$state$fitted_timegrid
      eta_T_mu <- eta_T_mu + x$mu$smooth.construct[[j]]$state$fitted_T
    }
  }
  # alpha
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- 
        drop(x$alpha$smooth.construct[[j]]$Xgrid %*% b)
      eta_timegrid_alpha <- eta_timegrid_alpha + 
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid
    }
  }
  # lambda
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <- 
        drop(x$lambda$smooth.construct[[j]]$Xgrid %*% b)
      eta_timegrid_lambda <- eta_timegrid_lambda + 
        x$lambda$smooth.construct[[j]]$state$fitted_timegrid
    }
  }
  
  # Eta predictors
  eta_timegrid_long <- rowSums(matrix(eta_timegrid_alpha*eta_timegrid_mu, 
                                      nrow = nsubj*nw, ncol = nmarker))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  eta_T_long <- rowSums(matrix(eta$alpha*eta_T_mu, nrow = nsubj,
                               ncol = nmarker))
  eta_T <- eta$lambda + eta$gamma + eta_T_long
  
  # Old logLikelihood and prior.
  sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
    (diag(nsubj)%x%t(gq_weights))%*% 
    exp(eta_timegrid)
  logLik_old <- drop(status %*% eta_T - sum_Lambda) +
    sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
              log = TRUE))
  
  # Fct for saving acceptance probability on the right scale
  transform_acceptprop <- function(x) {
    if(is.null(x)) return(0)
    if(is.na(x)) return(0)
    x <- exp(x)
    if(x < 0)
      x <- 0
    if(x > 1)
      x <- 1
    x
  }
  
  
  ## Process iterations
  if (burnin < 1) burnin <- 1
  if (burnin > n.iter) burnin <- floor(n.iter * 0.1)
  if (thin < 1) thin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))
  
  ## Samples.
  samps <- list()
  for (i in nx) {
    samps[[i]] <- list()
    for(j in names(x[[i]]$smooth.construct)) {
      samps[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), 
          ncol = length(x[[i]]$smooth.construct[[j]]$state$parameters)),
        "edf" = rep(NA, length = length(iterthin)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(samps[[i]][[j]]$samples) <- 
        names(x[[i]]$smooth.construct[[j]]$state$parameters)
    }
  }
  logLik.samps <- logPost.samps <- rep(NA, length = length(iterthin))

  nstep <- step
  step <- floor(n.iter / step)
  
  ptm <- proc.time()
  for(iter in 1:n.iter) {
    
    if (!is.null(prop_list)) {
      nx_iter <- 1
    }
    if(save <- iter %in% iterthin) {
      js <- which(iterthin == iter)
    }
    
    if(verbose_sampler) {
      cat("Iteration", iter, "\n")
    }
    #if(iter == 2) browser()
    
    for (i in nx) {
      
      if(!is.null(prop_list)) {
        j_iter <- 1
      }
      for (j in names(x[[i]]$smooth.construct)) {
        p_state <- propose_mjm(predictor = i,
                               x = x[[i]]$smooth.construct[[j]], y = y,
                               eta = eta, eta_timegrid = eta_timegrid,
                               eta_T = eta_T, eta_T_mu = eta_T_mu,
                               eta_timegrid_alpha = eta_timegrid_alpha,
                               eta_timegrid_mu = eta_timegrid_mu, 
                               eta_timegrid_long = eta_timegrid_long,
                               eta_timegrid_lambda = eta_timegrid_lambda,
                               survtime = survtime, logLik_old = logLik_old, 
                               nsubj = nsubj, gq_weights = gq_weights, 
                               status = status, nmarker = nmarker, 
                               nu = nu_sampler,
                               verbose_sampler = verbose_sampler, 
                               prop = if(!is.null(prop_list)) {
                                 prop_list[[iter]][[nx_iter]][[j_iter]]
                               } else NULL)
        
        # If accepted, set current state to proposed state
        accepted <- if(!is.null(prop_list)) {
          attr(prop_list[[iter]][[nx_iter]][[j_iter]], "acc")
        } else {if(is.na(p_state$xstate$alpha)){
          FALSE
        } else {
          log(runif(1)) <= p_state$xstate$alpha
        }}
        #cat(i, ": ", j, " - ", accepted, "\n")
        if (accepted) {
          
          # Update the etas
          switch(i, "lambda" = {
            eta_T <- p_state$etas$eta_T
            eta_timegrid <- p_state$etas$eta_timegrid 
            eta_timegrid_lambda <- p_state$etas$eta_timegrid_lambda
          }, "gamma" = {
            eta_T <- p_state$etas$eta_T
          }, "alpha" = {
            eta_T <- p_state$etas$eta_T
            eta_timegrid <- p_state$etas$eta_timegrid
            eta_timegrid_long <- p_state$etas$eta_timegrid_long
            eta_timegrid_alpha <- p_state$etas$eta_timegrid_alpha
          }, "mu" = {
            eta_T <- p_state$etas$eta_T
            eta_T_mu <- p_state$etas$eta_T_mu
            eta_timegrid <- p_state$etas$eta_timegrid
            eta_timegrid_long <- p_state$etas$eta_timegrid_long
            eta_timegrid_mu <- p_state$etas$eta_timegrid_mu
          })
          eta <- p_state$etas$eta
          
          # Update likelihood and state
          logLik_old <- p_state$logLik
          x[[i]]$smooth.construct[[j]]$state <- p_state$xstate
        } 
        #### SOLLTEN NICHT DIE VARIANZ-PARAMETER IMMER UPGEDATED WERDEN?
        #else {
        #   x[[i]]$smooth.construct[[j]]$state$parameters <- 
        #     bamlss::set.par(
        #       x[[i]]$smooth.construct[[j]]$state$parameters, 
        #       bamlss::get.par(x[[i]]$smooth.construct[[j]]$state$parameters, 
        #                       "tau2"), 
        #       "tau2")
        # }
        
        ## Save the samples and acceptance.
        if(save) {
          samps[[i]][[j]]$samples[js, ] <- 
            x[[i]]$smooth.construct[[j]]$state$parameters
          samps[[i]][[j]]$edf[js] <- x[[i]]$smooth.construct[[j]]$state$edf
          samps[[i]][[j]]$alpha[js] <- 
            transform_acceptprop(p_state$xstate$alpha)
          samps[[i]][[j]]$accepted[js] <- accepted
        }
        if (!is.null(prop_list)){
          j_iter <- j_iter + 1
        }
      }
      
      if(!is.null(prop_list)) {
        nx_iter <- nx_iter + 1
      }
    }
    if(save) {
      logLik.samps[js] <- logLik_old
      logPost.samps[js] <- as.numeric(logLik.samps[js] + 
                                        bamlss:::get.log.prior(x))
    }
  }
  
  for(i in names(samps)) {
    for(j in names(samps[[i]])) {
      samps[[i]][[j]] <- do.call("cbind", samps[[i]][[j]])
      cn <- if(j == "model.matrix") {
        paste(i, "p", j, colnames(samps[[i]][[j]]), sep = ".")
      } else {
        paste(i, "s", j, colnames(samps[[i]][[j]]), sep = ".")
      }
      colnames(samps[[i]][[j]]) <- cn
    }
    samps[[i]] <- do.call("cbind", samps[[i]])
  }
  samps$logLik <- logLik.samps
  samps$logPost <- logPost.samps
  samps <- do.call("cbind", samps)
  
  ## Compute DIC. #
  dev <- -2 * logLik.samps
  mpar <- apply(samps, 2, mean, na.rm = TRUE)
  names(mpar) <- colnames(samps)
  ll <- family$p2logLik(mpar)
  mdev <- -2 * ll
  pd <- mean(dev) - mdev
  DIC <- mdev + 2 * pd
  samps <- cbind(samps,
                 "DIC" = rep(DIC, length.out = nrow(samps)),
                 "pd" = rep(pd, length.out = nrow(samps))
  )
  if (is.null(prop_pred)) {
    samps[is.na(samps)] <- 0
  }
  
  
  return(as.mcmc(samps))
}

