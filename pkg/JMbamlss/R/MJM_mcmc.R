
# MJM_mcmc ----------------------------------------------------------------

MJM_mcmc <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
                     n.iter = 1200, burnin = 200, thin = 1, step = 20, 
                     nu_sampler = 1, ...)
{
  
  # Set starting values for the sampling
  if(!is.null(start)) {
    if(is.matrix(start)) {
      if(any(i <- grepl("Mean", colnames(start))))
        start <- start[, i]
      else stop("the starting values should be a vector not a matrix!")
    }
    x <- bamlss:::set.starting.values(x, start)
  }
  
  ## Names of parameters/predictors.
  nx <- names(x)
  nmarker <- attr(y, "nmarker")
  take_last <- attr(y, "take_last")
  survtime <- y[[1]][, "time"][take_last]
  nsubj <- length(survtime)
  
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
  
  eta_timegrid_long <- drop(
    t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
      (eta_timegrid_alpha*eta_timegrid_mu))
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  
  
  
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
  
  
  # UMBENENNEN: Warum wird das hier gebraucht? Woher kommt die exp()Fkt?
  foo <- function(x) {
    if(is.null(x)) return(0)
    if(is.na(x)) return(0)
    x <- exp(x)
    if(x < 0)
      x <- 0
    if(x > 1)
      x <- 1
    x
  }
  
  ## Integrals.
  # eeta <- exp(eta_timegrid)
  # int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + 
  #         apply(eeta[, 2:(sub - 1)], 1, sum))
  
  
  nstep <- step
  step <- floor(n.iter / step)
  
  ptm <- proc.time()
  for(iter in 1:n.iter) {
    if(save <- iter %in% iterthin)
      js <- which(iterthin == iter)
    
    #for(i in nx) {
      # if(i == "gamma") {
      #   eeta <- exp(eta_timegrid)
      #   int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + 
      #                      apply(eeta[, 2:(sub - 1)], 1, sum))
      # }
      
      # prop_fun <- get_jm_prop_fun(i, slice[i], nonlinear)
      # 3jump
    
    eta_T_long <- drop(
      t(rep(1, nmarker)) %x% diag(nsubj) %*% (eta$alpha*eta_T_mu))
    eta_T <- eta$lambda + eta$gamma + eta_T_long
    
    ## (1) update lambda.
    for(j in names(x[[i]]$smooth.construct)) {
      p.state <- propose_mjm_lambda(x = x[[i]]$smooth.construct[[j]], y = y,
                                    eta = eta, eta_timegrid = eta_timegrid,
                                    eta_timegrid_lambda = eta_timegrid_lambda,
                                    eta_T = eta_T, survtime = survtime, 
                                    nu = nu_sampler)
      eta_timegrid_lambda <- eta_timegrid_lambda -
        x$lambda$smooth.construct[[j]]$state$fitted_timegrid +
        state$fitted_timegrid
      eta_timegrid <- eta_timegrid_lambda +
        eta_timegrid_long
      eta$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[j]]$state) +
        fitted(state)
      x$lambda$smooth.construct[[j]]$state <- state
    }
    
    
    
      for(sj in names(x[[i]]$smooth.construct)) {
        
        p.state <- prop_fun(x[[i]]$smooth.construct[[sj]],
                            y, eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                            eta_timegrid_dmu, eta_timegrid_dalpha,
                            width, sub, nu, status, id = i, int0, nobs,
                            dx = if(dalpha & (i == "mu")) x[["dmu"]]$smooth.construct[[sj]] else NULL,
                            xsmalpha = if(nonlinear & (i == "mu")) x$alpha$smooth.construct else NULL, 
                            knots = if(nonlinear & (i == "mu")) knots else NULL, tp = tp, fac = fac, ...)
        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha
        
        #checks# if(i == "mu") cat("\n", sj, round(exp(p.state$alpha), 2))
        if(i %in% fixed)
          accepted <- FALSE
        
        if(accepted) {
          if(i %in% c("lambda", "mu", "alpha", "dalpha")) {
            if(i == "lambda")
              eta_timegrid_lambda <- eta_timegrid_lambda - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            if(i == "mu") {
              # cat("\n iteration: ", iter, ", predictor: ", i, ", term: ", sj)
              eta_timegrid_mu <- eta_timegrid_mu - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
              if(nonlinear){
                for (j in names(x$alpha$smooth.construct)){
                  if(j != "model.matrix"){   # only smooth.constructs need to be updated
                    g_a <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
                    Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
                    Xalpha <- constrain(x$alpha$smooth.construct[[j]], Xalpha)
                    if(fac)
                      Xalpha <- Xalpha * x$alpha$smooth.construct[[j]]$by_timegrid
                    if(tp)  
                      Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
                    alpha_state <- matrix(Xalpha %*% g_a, nrow = nrow(eta_timegrid), ncol = ncol(eta_timegrid), byrow = TRUE)
                    eta_timegrid_alpha <- eta_timegrid_alpha - x$alpha$smooth.construct[[j]]$state$fitted_timegrid + alpha_state
                    x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- alpha_state
                    eta$alpha <- eta$alpha - fitted(x$alpha$smooth.construct[[j]]$state) + alpha_state[, ncol(eta_timegrid)]
                    x$alpha$smooth.construct[[j]]$state$fitted.values <- alpha_state[, ncol(eta_timegrid)]
                  } 
                }
              }
              if(dalpha & (sj %in% names(x[["dmu"]]$smooth.construct))) {
                p.state.dmu <- update_jm_dmu(x[["dmu"]]$smooth.construct[[sj]], x[["mu"]]$smooth.construct[[sj]])
                eta_timegrid_dmu <- eta_timegrid_dmu - x[["dmu"]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state.dmu$fitted_timegrid
                eta[["dmu"]] <- eta[["dmu"]] - fitted(x[["dmu"]]$smooth.construct[[sj]]$state) + fitted(p.state.dmu)
                x[["dmu"]]$smooth.construct[[sj]]$state <- p.state.dmu
              }
            }
            if(i == "alpha"){
              eta_timegrid_alpha <- eta_timegrid_alpha - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            }
            if(i == "dalpha")
              eta_timegrid_dalpha <- eta_timegrid_dalpha - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            if(nonlinear){
              eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
            } else {
              eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
            }
          }
          eta[[i]] <- eta[[i]] - fitted(x[[i]]$smooth.construct[[sj]]$state) + fitted(p.state)
          if(i == "alpha"){
            if(myplots & iter == round(iter, -1)){
              par(mfrow=c(2,2))
              plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)))
              plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha, main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)))
              plot(eta_timegrid_mu, p.state$fitted_timegrid, 
                   main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2), ", term: ", sj))
              abline(h = 0, col = "red")
              abline(v = median(y[, 3]), col = "red")
              matplot(t(times), t(eta_timegrid_alpha), main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)), type = "l")
              Sys.sleep(2)
            }
          }
          if(i == "mu"){
            if(myplots & iter == round(iter, -1)){
              par(mfrow=c(2,2))
              matplot(t(times),t(p.state$fitted_timegrid), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), " ", sj))
              matplot(t(times), t(eta_timegrid_mu), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), " ", sj))
              plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), "alpha"))
              plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha, main=paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), "alpha"))
              Sys.sleep(1)
            }
          }
          x[[i]]$smooth.construct[[sj]]$state <- p.state 
        }
        
        ## Save the samples and acceptance.
        if(save) {
          samps[[i]][[sj]]$samples[js, ] <- x[[i]]$smooth.construct[[sj]]$state$parameters
          samps[[i]][[sj]]$edf[js] <- x[[i]]$smooth.construct[[sj]]$state$edf
          samps[[i]][[sj]]$alpha[js] <- foo(p.state$alpha)
          samps[[i]][[sj]]$accepted[js] <- accepted
          if(dalpha & (i == "mu") & (sj %in% names(x[["dmu"]]$smooth.construct))) {
            samps[["dmu"]][[sj]]$samples[js, ] <- x[["dmu"]]$smooth.construct[[sj]]$state$parameters
            samps[["dmu"]][[sj]]$edf[js] <- x[["dmu"]]$smooth.construct[[sj]]$state$edf
            samps[["dmu"]][[sj]]$alpha[js] <- foo(p.state$alpha)
            samps[["dmu"]][[sj]]$accepted[js] <- accepted
          }
        }
      }
    #}
    
    if(save) {
      logLik.samps[js] <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
        exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
      logPost.samps[js] <- as.numeric(logLik.samps[js] + get.log.prior(x))
    }
    
    if(verbose) barfun(ptm, n.iter, iter, step, nstep)
  }
  
  if(verbose) cat("\n")
  
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
  
  return(as.mcmc(samps))
}

