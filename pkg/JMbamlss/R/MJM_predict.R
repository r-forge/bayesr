## Prediction.
MJM_predict <- function(object, newdata,
                        type = c("link", "parameter", "probabilities", "cumhaz",
                                 "loglik"),
                       dt, steps, id,
                       FUN = function(x) { mean(x, na.rm = TRUE) },
                       subdivisions = 7, cores = NULL,
                       chunks = 1, verbose = FALSE,  ...)
{
  # if(attr(object$y[[1]], "nonlinear") & 
  #    any(type %in% c("probabilities", "cumhaz"))){
  #   stop("Computation of cumulative hazard and survival probabilities", 
  #        "are currently under construction for associations which are",
  #        "nonlinear in mu."))
  # }
  if(missing(dt)) dt <- 0
  if(missing(steps)) steps <- 1
  if(missing(id)) i <- NULL else i <- id
  idvar <- attr(object$y, "idvar")
  marker_name <- attr(object$y, "marker_name")
  nmarker <- attr(object$y, "nmarker")
  marker_levels <- levels(attr(object$y, "marker"))
  timevar_mu <- attr(object$y, "timevar")["mu"]
  
  if(length(type) > 1)
    type <- type[1]
  type <- match.arg(type)
  loglik <- type == "loglik"
  if(loglik) {
    type <- "cumhaz"
    dt <- 0
  }
  
  if(type == "probabilities"){
    if(!is.null(newdata)){
      warning("The provided newdata will be ignored for the prediction",
              "of conditional survival probabilities.")
      newdata <- NULL
    }
    if(dt == 0){
      stop("Please specify a time window for the prediction of conditional", 
           "survival probabilities.")
    }
  }
  
  if(type == "cumhaz" & !is.null(newdata) & dt > 0){
    warning("The provided newdata will be ignored for the prediction of", 
            "conditional survival t + dt.")
    newdata <- NULL
  }
  
  # calculate_mu <- FALSE
  # if(
  #   # only if nonlinear effect
  #   attr(object$y[[1]], "nonlinear") &
  #   # only if alpha amongst predicted effects
  #   (("alpha" %in% list(...)$model) | is.null(list(...)$model)) & 
  #   # only if mu is not explicitely specified in newdata
  #   (is.null(newdata)|is.null(newdata$mu))                        
  # ){
  #   calculate_mu <- TRUE
  # }
  
  if(is.null(newdata)){
    newdata <- model.frame(object) 
  }
  
  # if(calculate_mu){
  #   object$family$predict <- NULL
  #   long_timevar <- attr(object$y[[1]], "timevar")["mu"]
  #   surv_timevar <- attr(object$y[[1]], "timevar")["lambda"]
  #   tempdata <- newdata
  #   tempdata[, long_timevar] <- tempdata[, surv_timevar]
  #   mu <- predict.bamlss(object, model = "mu", newdata = tempdata, 
  #                        type = type, cores = cores, chunks = chunks,
  #                        verbose = verbose)
  #   newdata$mu <- mu
  # }
  
  if(!(type %in% c("probabilities", "cumhaz"))) {
    
    object$family$predict <- NULL
    
    # Check for PCRE-term
    pcre_mu <- lapply(object$x$mu$smooth.construct, function (x) {
      any(grepl("pcre.random.effect", class(x), fixed = TRUE))
    })
    for (j in which(unlist(pcre_mu))) {
      
      # Compute the evaluated fpc basis functions for the newdata
      attr(object$x$mu$smooth.construct[[j]], "which_marker") <- 
        newdata[, attr(object$y, "marker_name")]
      
      # Check and attach missing FPCbasis (set to 0)
      add_fpc <- which(!object$x$mu$smooth.construct[[j]]$term %in% 
                         colnames(newdata))
      if (length(add_fpc) > 0) {
        for (add in add_fpc) {
          newdata[[object$x$mu$smooth.construct[[j]]$term[add]]] <- 0
        }
      }
      
    }
    
    return(bamlss:::predict.bamlss(object, newdata = newdata, type = type,
                                   FUN = FUN, cores = cores, chunks = chunks, 
                                   verbose = verbose, ...))
  }
  
  if(object$family$family != "mjm")
    stop("object must be a mjm-type model!")
  
  # dalpha <- has_pterms(object$x$dalpha$terms) | 
  #   (length(object$x$dalpha$smooth.construct) > 0)
  
  timevar <- attr(object$y, "timevar")
  tmax_model <- max(newdata[,timevar["lambda"]])
  
  if(!is.null(i)){
    if(!is.character(i))
      i <- levels(newdata[[idvar]])[i]
    newdata <- subset(newdata, newdata[[idvar]] %in% i)
  }
  
  tmax_pred <- max(newdata[,timevar["mu"]]) + dt
  
  if(tmax_pred > tmax_model){
    warning("Predictions should not be made beyond the modelled time range. 
            Please adjust the time window dt accordingly.")
  }
  
  
  
  ## Create the time grid. 
  # Gaussian Quadrature
  stopifnot(requireNamespace("statmod"))
  gq <- statmod::gauss.quad(subdivisions)
  
  grid <- function(lower, upper) {
    (upper - lower) / 2 * gq$nodes + (upper + lower) / 2
  }
  
  jm_probs <- function(data) {
    
    if(dt == 0){
      take <- !duplicated(data[, c(timevar["lambda"], idvar)])
      dsurv <- subset(data, take)
      timegrid <- lapply(dsurv[[timevar["lambda"]]], 
                         function(x){grid(0, x)})
    } else {
      take <- !duplicated(data[, c(timevar["lambda"], idvar)], fromLast = TRUE)
      dsurv <- subset(data, take)
      timegrid <- lapply(dsurv[[timevar["mu"]]], 
                         function(x){grid(x, x+dt)})
    }
    nobs <- nrow(dsurv)
    gdim <- c(length(timegrid), length(timegrid[[1]]))
    
    # long data grid for multiple markers
    # dsurv_long <- dsurv[rep(seq_len(nrow(dsurv)), nmarker), ]
    # dsurv_long[[marker_name]] <- rep(marker_levels, each = nrow(dsurv))
    # timegrid_long <- rep(timegrid, nmarker)
    # width <- rep(NA, nobs)
    # for(i in 1:nobs)
    #   width[i] <- timegrid[[i]][2] - timegrid[[i]][1]
    
    pred.setup <- bamlss:::predict.bamlss(object, data, type = "link",
                                          get.bamlss.predict.setup = TRUE, ...)
    
    enames <- pred.setup$enames
    
    pred_gamma <- with(pred.setup, 
                       bamlss:::.predict.bamlss(
                         "gamma", object$x$gamma, samps, enames$gamma, 
                         intercept, nsamps, dsurv))
    
    pred_lambda <- with(pred.setup, 
                        .predict.bamlss.mjm.td(
                          id = "lambda", 
                          x = object$x$lambda$smooth.construct, samps = samps,
                          enames = enames$lambda, intercept = intercept, 
                          nsamps = nsamps, newdata = dsurv, 
                          # yname and timevar needed, take and y as well
                          yname = timevar["lambda"], grid = timegrid, 
                          formula = bamlss:::drop.terms.bamlss(
                            object$x$lambda$terms, sterms = FALSE, 
                            keep.response = FALSE), idvar = idvar, 
                          timevar_mu = timevar_mu, nmarker = nmarker))
    # type = 2 bedeutete hier nur, dass param_time_transform2 verwendet wird

    #!# Todo: Fix prediction for nonlinear
    # Adapt grid -> longer, for each marker
    # Adapt timevar to timevar_mu
    pred_mu <- with(pred.setup, 
                    .predict.bamlss.mjm.td(
                      id = "mu", x = object$x$mu$smooth.construct, 
                      samps = samps, enames = enames$mu, intercept = intercept,
                      nsamps = nsamps, newdata = dsurv, # dsurv_long, 
                      yname = timevar["lambda"], 
                      grid = timegrid, 
                      formula = bamlss:::drop.terms.bamlss(
                        object$x$mu$terms, sterms = FALSE, 
                        keep.response = FALSE), idvar = NULL, 
                      timevar_mu = timevar_mu, nmarker = nmarker))
    
    pred_alpha <- with(pred.setup, 
                       .predict.bamlss.mjm.td(
                         id = "alpha", object$x$alpha$smooth.construct, 
                         samps = samps, enames = enames$alpha, 
                         intercept = intercept, nsamps = nsamps,
                         newdata = dsurv,  yname = timevar["lambda"], 
                         grid = timegrid_long,
                         formula = bamlss:::drop.terms.bamlss(
                           object$x$alpha$terms, sterms = FALSE,
                           keep.response = FALSE), idvar = NULL, 
                         timevar_mu = timevar_mu, nmarker = nmarker))
    
    # if(dalpha) {
    #   pred_dalpha <- with(pred.setup, 
    #                       .predict.bamlss.surv.td(
    #                         "dalpha", object$x$dalpha$smooth.construct, samps,
    #                         enames$dalpha, intercept, nsamps, dsurv,
    #                         timevar["lambda"], timegrid, 
    #                         drop.terms.bamlss(object$x$dalpha$terms,
    #                                           sterms = FALSE, 
    #                                           keep.response = FALSE)))
    #   
    #   pred_dmu <- with(pred.setup, 
    #                    .predict.bamlss.surv.td(
    #                      "dmu", object$x$dmu$smooth.construct, samps, 
    #                      enames$dmu, intercept, nsamps, dsurv, timevar["mu"], 
    #                      timegrid, 
    #                      drop.terms.bamlss(object$x$dmu$terms, sterms = FALSE,
    #                                        keep.response = FALSE),
    #                      derivMat = TRUE))
    # }
    # if(attr(object$y[[1]], "nonlinear")){
    #   eta_timegrid <- if(dalpha) {
    #     pred_lambda + pred_alpha + pred_dalpha * pred_dmu
    #   } else {
    #     pred_lambda + pred_alpha
    #   }
    # } else {
      eta_timegrid <- #if(dalpha) {
      #   pred_lambda + pred_alpha * pred_mu + pred_dalpha * pred_dmu
      # } else {
        pred_lambda + pred_alpha * pred_mu
      #}
    # }
    
    if(loglik) {
      eta_gamma <- bamlss:::predict.bamlss(object, data[take, , drop = FALSE], 
                                           model = "gamma")
      eta_mu <- bamlss:::predict.bamlss(object, data, model = "mu")
      eta_sigma <- bamlss:::predict.bamlss(object, data, model = "sigma")
      mf <- model.frame(object, data = data)
      y <- mf[, grep("Surv", names(mf), fixed = TRUE)]
      eta_timegrid <- matrix(eta_timegrid, nrow = gdim[1], ncol = gdim[2], 
                             byrow = TRUE)
      eeta <- exp(eta_timegrid)
      # ANPASSEN
      int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + 
                        apply(eeta[, 2:(subdivisions - 1)], 1, sum))
      ll <- sum((eta_timegrid[, ncol(eta_timegrid)] + 
                   eta_gamma) * y[take, "status"], na.rm = TRUE) -
        exp(eta_gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta_mu,
                                           sd = exp(eta_sigma), log = TRUE))
      return(drop(ll))
    }
    
    if(dt == 0){
      probs <- NULL
      for(i in 1:ncol(eta_timegrid)) {
        eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2],
                      byrow = TRUE)
        eeta <- exp(eta)
        # ANPASSEN
        int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + 
                          apply(eeta[, 2:(subdivisions - 1), drop = FALSE], 1, 
                                sum))
        probs <- if(type == "probabilities") {
          cbind(probs, exp(-1 * exp(pred_gamma[, i]) * int))
        } else {
          cbind(probs, exp(pred_gamma[, i]) * int)
        }
      }
      
      if(!is.null(FUN)) {
        if(is.matrix(probs)) {
          if(ncol(probs) > 1)
            probs <- apply(probs, 1, FUN)
          probs <- t(probs) 
        } 
      }
      
    } else {
      lprobs <- lapply(1:steps, function(x){
        probs <- NULL
        sub <- round(subdivisions/steps * x)
        for(i in 1:ncol(eta_timegrid)) {
          eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2], 
                        byrow = TRUE)
          eeta <- exp(eta)
          # ANPASSEN
          int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + 
                            apply(eeta[, 2:(sub - 1), drop = FALSE], 1, sum))
          probs <- if(type == "probabilities") {
            cbind(probs, exp(-1 * exp(pred_gamma[, i]) * int))
          } else {
            cbind(probs, exp(pred_gamma[, i]) * int)
          }
        }
        
        if(!is.null(FUN)) {
          if(is.matrix(probs)) {
            if(ncol(probs) > 1)
              probs <- apply(probs, 1, FUN)
            probs <- t(probs) 
          } 
        }
        probs
      }) 
      
      probs <- lprobs
      names(probs) <- paste("Time after last longitudinal observation:", 
                            (1:steps)*dt/steps) 
    }
    
    
    return(probs)
  }
  #debug(jm_probs)
  ia <- interactive()
  
  if(is.null(cores)) {
    if(chunks < 2) {
      probs <- jm_probs(newdata)
    } else {
      id <- sort(rep(1:chunks, length.out = nrow(newdata)))
      newdata <- split(newdata, id)
      chunks <- length(newdata)
      probs <- NULL
      for(i in 1:chunks) {
        if(verbose) {
          cat(if(ia) "\r" else "\n")
          cat("predicting chunk", i, "of", chunks, "...")
          if(.Platform$OS.type != "unix" & ia) flush.console()
        }
        if(i < 2) {
          probs <- jm_probs(newdata[[i]])
        } else {
          if(is.null(dim(probs))) {
            probs <- c(probs, jm_probs(newdata[[i]]))
          } else {
            probs <- rbind(probs, jm_probs(newdata[[i]]))
          }
        }
      }
      if(verbose) cat("\n")
    }
  } else {
    parallel_fun <- function(i) {
      if(chunks < 2) {
        pr <- jm_probs(newdata[[i]])
      } else {
        idc <- sort(rep(1:chunks, length.out = nrow(newdata[[i]])))
        nd <- split(newdata[[i]], idc)
        chunks <- length(nd)
        pr <- NULL
        for(j in 1:chunks) {
          if(j < 2) {
            pr <- jm_probs(nd[[j]])
          } else {
            if(is.null(dim(pr))) {
              pr <- c(pr, jm_probs(nd[[j]]))
            } else {
              pr <- rbind(pr, jm_probs(nd[[j]]))
            }
          }
        }
      }
      return(pr)
    }
    
    id <- sort(rep(1:cores, length.out = nrow(newdata)))
    newdata <- split(newdata, id)
    cores <- length(newdata)
    probs <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
    
    probs <- if(is.matrix(probs[[1]])) {
      do.call("rbind", probs)
    } else {
      do.call("c", probs)
    }
  }
  
  return(probs)
}


# Vielleicht muss man hier eine eigene Version schreiben?
.predict.bamlss.mjm.td <- function(id, x, samps, enames, intercept, nsamps, 
                                   newdata, yname, grid, formula, idvar, 
                                   timevar_mu, nmarker)
{
 
  # id ist "lambda", "mu" etc, also vielleicht umbenennen
  snames <- colnames(samps)
  enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
  has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."), 
                             snames, fixed = TRUE))
  
  # Warum sollte man hier noch einen zusätzlichen Intercept brauchen?
  if(intercept & has_intercept)
    enames <- c("p.(Intercept)", enames)
  if (!has_intercept) {
    enames <- enames[-grep("p.(Intercept)", enames, fixed = TRUE)]
  }
  enames <- unique(enames)
  ec <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][1:2], collapse = "")
  })
  enames2 <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][-c(1:2)], collapse = "")
  })
  
  eta <- 0
  p_components <- grep("p.", ec)
  if(length(p_components)) {
    intcpt <- ifelse(has_intercept, "1", "-1")
    f <- as.formula(paste("~", intcpt, "+",
                          paste(enames2[p_components], collapse = "+")))
    X <- param_time_transform_mjm(x = list(), formula = f, data = newdata, 
                                  grid = grid, yname = yname,
                                  timevar = yname, take = NULL, 
                                  idvar = idvar, y, timevar2 = 
                                    if (i == lambda) {timevar_mu} 
                                  else {NULL})$Xgrid
    
    sn <- snames[bamlss:::grep2(paste(id, "p", 
                                      enames2[p_components], sep = "."), 
                                snames, fixed = TRUE)]
    if(!length(sn))
      sn <- snames[bamlss:::grep2(paste(id, "p.model.matrix", 
                                        enames2[p_components], sep = "."),
                                  snames, fixed = TRUE)]
    eta <- eta + bamlss:::fitted_matrix(X, samps[, sn, drop = FALSE])
  }
  if(length(i <- grep("s.", ec))) {
    y2 <- newdata[, yname, drop = FALSE]
    for(j in enames2[i]) {
      for(jj in grep(j, names(x), fixed = TRUE, value = TRUE)) {
        
        if(!inherits(x[[jj]], "no.mgcv") & !inherits(x[[jj]], "special")) {
          if(inherits(x[[j]], 
                      "pcre.random.effect")) {
            browser()
            x[[jj]]$term <- x[[jj]]$term[-length(x[[jj]]$term)] 
            X <- sm_time_transform_mjm_pcre(
              x = x[[jj]], data = newdata, 
              grid = grid, yname = yname, timevar = timevar_mu,
              take = NULL, nmarker = nmarker)$Xgrid
          } else {
            X <- sm_time_transform_mjm(
              x = x[[jj]], data = newdata, 
              grid = grid, yname = yname, 
              timevar = if (id == "mu") timevar_mu else yname,
              take = NULL, y = y2)$Xgrid
          }
          
          # ändern zu sm_time_transform_mjm
          #take_last <- rep(TRUE, nrow(newdata))
          
          
          # X <- sm_Xtimegrid(x[[jj]], newdata, grid, yname,
          #                   derivMat = derivMat)
          sn <- snames[bamlss:::grep2(paste(id, "s", jj, sep = "."), snames, 
                             fixed = TRUE)]
          random <- if(!is.null(x[[jj]]$margin)) {
            any(sapply(x[[jj]]$margin, function(z) {
              inherits(z, "random.effect") 
              }))
          } else inherits(x[[jj]], "random.effect")
          ok <- if(random) {
            if(ncol(X) == length(samps[, sn, drop = FALSE])) TRUE else FALSE
          } else TRUE
          if(ok) {
            eta <- if(is.null(eta)) {
              bamlss:::fitted_matrix(X, samps[, sn, drop = FALSE])
            } else {
              eta + bamlss:::fitted_matrix(X, samps[, sn, drop = FALSE])
            }
          }
        } else {
          stop("no predictions for special terms available yet!")
        }
      }
    }
  }
  
  eta
}
