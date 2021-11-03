
if(!exists("d")) {
  library("refund")
  library("bamlss")
  library("MFPCA")
  
  source("simMultiJM.R")
  source("eval_mfun.R")
  
  d <- simMultiJM(nsub = 50)

  marker_dat <- split(d, d$marker)
  marker_dat <- lapply(marker_dat, function (mark) {
    mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
    mark
  })
  m_irregFunData <- lapply(marker_dat, function (mark) {
    mark <- mark[order(mark$obstime), ]
    irregFunData(argvals = split(mark$obstime, mark$id), 
                 X = split(mark$res, mark$id))
  })
  FPCA <- lapply(m_irregFunData, function(mark) {
    PACE(mark)
  })
  mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
  uniExpansions <- lapply(FPCA, function (mark) {
    list(type = "given", functions = mark$functions)
  })
  MFPCA <- MFPCA(mFData = mFData, M = 2, uniExpansions = uniExpansions)
 
}

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
    "sampler" = NULL,
    "predict" = NULL
  )
  
  class(rval) <- "family.bamlss"
  rval
}

## Smooth time transformer function.
sm_time_transform_mjm <- function(x, data, grid, yname, timevar, take, y)
{

  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- XT <- NULL
  for(j in x$term) {
    if((j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
      XT <- if(is.null(XT)) df else cbind(XT, df)
    }
  }
  if(!is.null(X)) {
    colnames(X) <- colnames(XT) <- x$term[!(x$term %in% c(yname, timevar))]
  }
    
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  XT <- if(is.null(XT)) data.frame(y[, 1]) else cbind(X, y[, 1])
  colnames(X)[ncol(X)] <- colnames(XT)[ncol(XT)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    XT <- cbind(XT, y[, 1])
    colnames(X)[ncol(X)] <- colnames(XT)[ncol(XT)] <- timevar
  }
  if(x$by != "NA" & x$by != yname) {
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))
    XT[[x$by]] <- data[[x$by]]
    if (nrow(XT) != nrow(data)) {
      stop("XT dimensions do not coincide with 'by' dimensions for ", x$label)
    }
  }
    
  x$Xgrid <- PredictMat(x, X)
  
  # XT necessary for calculation of score and hessian
  # Matrix of evaluations at the vector of survival times
  x$XT <- PredictMat(x, XT)
  
  x
}

## PCRE transformer.
sm_time_transform_mjm_pcre <- function(x, data, grid, yname, timevar, take,
                                       nmarker) {
  
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- NULL
  for(j in x$term) {
    if((j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
    }
  }
  if(!is.null(X))
    colnames(X) <- x$term[!(x$term %in% c(yname, timevar))]
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
    data[[timevar]] <- data[[yname]]
  }
  if(x$by != "NA" & x$by != yname)
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))

  class(x) <- "pcre2.random.effect"
  x$term <- c(x$term, timevar)
  x$timevar <- timevar
  x$Xgrid <- PredictMat(x, X, n = nmarker*nrow(X))
  x$XT <- PredictMat(x, data, n = nmarker*nrow(data))
  
  x
}

Predict.matrix.pcre2.random.effect <- function(object, data)
{
  if(is.null(object$xt$mfpc))
    stop("need mfpa object!")
  X <- eval_mfpc(object$xt$mfpc, data[[object$timevar]])
  if(ncol(X) != (length(object$term) - 2))
    stop("check M argument in MFPCA()!")
  X <- data.frame(data[[object$term[1]]], X)
  colnames(X) <- object$term[-length(object$term)]
  
  # Muss man dann dieses X vielleicht noch als data argument in die
  # Predict.matrix.pcre.random.effect()
  # übergeben?
  object$term <- object$term[-length(object$term)]
  X <- refund:::Predict.matrix.pcre.random.effect(object, X)
  return(X)
}

## Linear design transformer.
param_time_transform_mjm <- function(x, formula, data, grid, yname, timevar, 
                                     take, idvar, y)
{

  X <- Xn <- tvar <- NULL
  # For time-varying covariates in lambda predictor (idvar is not NULL)
  if (!is.null(idvar)) {
    id <- data[[idvar]]
    for (j in names(data)) {
      if ((!grepl("Surv(", j, fixed = TRUE) & 
           !grepl("Surv2(", j, fixed = TRUE)) & (j != yname) & (j != timevar)) {
        
        # split data per subject
        idata <- split(data[[j]], id)
        # check if timevarying variable
        temp <- lapply(1:length(idata), function(i){
          length(unique(idata[[i]])) > 1
          })
        if(any(unlist(temp))){
          tvar <- c(tvar, j)
          # extract unique time-varying values
          values <- lapply(1:length(idata), function(i){unique(idata[[i]])})
          # extract break points
          breaks <- lapply(1:length(idata), function(i){
            split(data[[timevar2]], id)[[i]][c(TRUE, diff(idata[[i]]) != 0)]})
          # transfer break points to evaluation grid
          igrid <- lapply(1:length(idata), function(i){
            if(length(breaks[[i]]) > 1){
              g <- cut(grid[[i]], breaks[[i]], labels=FALSE,
                       include.lowest = TRUE)
              g[is.na(g)] <- max(g, na.rm=TRUE) + 1
              g
            } else {
              rep(1, length(grid[[i]]))
            }})
          # evaluate variable on that grid
          evalgrid <- lapply(1:length(idata), function(i){
            values[[i]][igrid[[i]]]})
          df <- data.frame(unlist(evalgrid))
          names(df) <- j
          X <- if (is.null(X)) 
            df
          else cbind(X, df)
          Xn <- c(Xn, j)
        }
      }
    }
  }
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  
  for(j in names(data)) {
    if((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", j, fixed = TRUE)) &
       (j != yname) & (j != timevar) & !(j %in% tvar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
      Xn <- c(Xn, j)
    }
  }
  if(!is.null(X)) {
    colnames(X) <- Xn
  }
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }
  
  x$Xgrid <- model.matrix(formula, data = X)
  x$XT <- model.matrix(formula, data = data)

  x
}


## Compute all necessary matrices.
MJM_transform <- function(object, subdivisions = 7, timevar = NULL, ...)
{
  stopifnot(requireNamespace("statmod"))

  gq <- statmod::gauss.quad(subdivisions, ...)

  ## Get idvar.
  class_mu <- sapply(object$x$mu$smooth.construct, class)
  class_mu <- sapply(class_mu, function(x) x[1])
  if(!any(class_mu == "pcre.random.effect"))
    stop("need pcre smooth!")
  j <- which(class_mu == "pcre.random.effect")
  smj <- object$x$mu$smooth.construct[[j]]
  idvar <- smj$term[1]

  ## Setup integration.
  create_grid <- function(time) {
    time / 2 * gq$nodes + time / 2
  }
  y2_l <- cbind(object$y[[1]][, "time"], object$model.frame[[idvar]])
  colnames(y2_l) <- c("time", idvar)

  
  take <- !duplicated(y2_l[, 1:2])
  take_last <- !duplicated(y2_l, fromLast = TRUE)
  nsubj <- length(unique(y2_l[, idvar]))
  y2 <- y2_l[take, , drop = FALSE]
  grid <- lapply(y2[, "time"], create_grid)
  

  nmarker <- FALSE
  marker <- rep(1, nrow(y2_l))
  if(!is.null(object$model.frame$marker)) {
    y2_l <- cbind(y2_l, "marker" = as.factor(object$model.frame$marker))
    nmarker <- length(levels(as.factor(object$model.frame$marker)))
    marker <- object$model.frame$marker
  }

  take_l <- !duplicated(y2_l)
  take_last_l <- !duplicated(y2_l, fromLast = TRUE)

  y2_l <- y2_l[take_l, , drop = FALSE]
  grid_l <- lapply(y2_l[, "time"], create_grid)

  ## Save information for optimizer in attributes of y
  attr(object$y, "gq_weights") <- gq$weights
  attr(object$y, "status") <- object$y[[1]][, "status"][take_last]
  attr(object$y, "take_last") <- take_last
  attr(object$y, "take_last_l") <- take_last_l
  attr(object$y, "nsubj") <- nsubj
  attr(object$y, "nmarker") <- nmarker
  attr(object$y, "marker") <- marker
  

  yname <- all.names(object$x$lambda$formula[2])[2]
  timevar_mu <- timevar
  if(is.null(timevar_mu))
    stop("the time variable is not specified, needed for mu!")
  timevar <- yname

  # design.construct als Funktion eingebaut, damit sich die eta-Vektoren
  # der Survival-Prädiktoren in der Länge ändern
  ## Recompute design matrixes for lambda, gamma, alpha.
  for(j in c("lambda", "gamma", "alpha")) {
    object$x[[j]] <- bamlss:::design.construct(object$terms, 
      data = object$model.frame[if (j == "alpha") take_l else take, , 
      drop = FALSE], model.matrix = TRUE, smooth.construct = TRUE, model = j, 
      scale.x = FALSE)[[j]]
  }
  
  ## The basic setup.
  if(is.null(attr(object$x, "bamlss.engine.setup")))
    object$x <- bamlss.engine.setup(object$x, ...)

  ## Remove intercept from lambda.
  if(!is.null(object$x$lambda$smooth.construct$model.matrix)) {
    cn <- colnames(object$x$lambda$smooth.construct$model.matrix$X)
    if("(Intercept)" %in% cn)
      object$x$lambda$smooth.construct$model.matrix$X <- 
        object$x$lambda$smooth.construct$model.matrix$X[, cn != "(Intercept)",
                                                        drop = FALSE]
    if(ncol(x$lambda$smooth.construct$model.matrix$X) < 1) {
      object$x$lambda$smooth.construct$model.matrix <- NULL
      object$x$lambda$terms <- drop.terms.bamlss(x$lambda$terms, pterms = FALSE,
                                                 keep.intercept = FALSE)
    } else {
      object$x$lambda$smooth.construct$model.matrix$term <- 
        gsub("(Intercept)+", "",
        object$x$lambda$smooth.construct$model.matrix$term, fixed = TRUE)
      object$x$lambda$smooth.construct$model.matrix$state$parameters <- 
        object$x$lambda$smooth.construct$model.matrix$state$parameters[-1]
      object$x$lambda$terms <- drop.terms.bamlss(x$lambda$terms,
        pterms = TRUE, sterms = TRUE, keep.intercept = FALSE)
    }
  }

  ## Compute new design matrices for integration for smooth terms.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct)) {
      for(j in names(object$x[[i]]$smooth.construct)) {
        if(j != "model.matrix") {
          xterm <- object$x[[i]]$smooth.construct[[j]]$term
          by <- if(object$x[[i]]$smooth.construct[[j]]$by != "NA") {
            object$x[[i]]$smooth.construct[[j]]$by
          } else NULL
          if(inherits(object$x[[i]]$smooth.construct[[j]], 
                      "pcre.random.effect")) {
            object$x[[i]]$smooth.construct[[j]] <- 
              sm_time_transform_mjm_pcre(object$x[[i]]$smooth.construct[[j]],
                object$model.frame[, unique(c(xterm, yname, by, timevar_mu,
                  idvar)), drop = FALSE], 
                grid, yname, timevar_mu, take_last, 
                nmarker = nmarker)
          } else {
            object$x[[i]]$smooth.construct[[j]] <- sm_time_transform_mjm(
              x = object$x[[i]]$smooth.construct[[j]],
              data = object$model.frame[, unique(c(xterm, yname, by,
                if(i == "mu") timevar_mu else timevar, idvar)), drop = FALSE],
              grid = if(i == "lambda") grid else grid_l,
              yname = yname, 
              timevar = if(i == "mu") timevar_mu else timevar,
              take = if(i == "lambda") take_last else take_last_l,
              y = if (i == "lambda") y2 else y2_l)
          }
        }
      }
    }
  }

  ## Now linear part.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct$model.matrix)) {
      object$x[[i]]$smooth.construct$model.matrix <- param_time_transform_mjm(
        object$x[[i]]$smooth.construct$model.matrix,
        bamlss:::drop.terms.bamlss(object$x[[i]]$terms, sterms = FALSE,
          keep.response = FALSE), object$model.frame,
        if(i == "lambda") grid else grid_l, yname, 
        if(i != "mu") timevar else timevar_mu,
        if(i == "lambda") take_last else take_last_l,
        idvar = if (i == "lambda") idvar else NULL)
    }
  }

  ## Update prior/grad/hess functions
  for(j in names(object$x)) {
    for(sj in names(object$x[[j]]$smooth.construct)) {
      priors <- bamlss:::make.prior(object$x[[j]]$smooth.construct[[sj]])
      object$x[[j]]$smooth.construct[[sj]]$prior <- priors$prior
      object$x[[j]]$smooth.construct[[sj]]$grad <- priors$grad
      object$x[[j]]$smooth.construct[[sj]]$hess <- priors$hess
    }
  }

  return(object)
}


opt_MJM <- function(x, y, start = NULL, eps = 0.0001, maxit = 600, nu = 0.1, ...)
{

  if(!is.null(start))
    x <- bamlss:::set.starting.values(x, start)
  
  nmarker <- attr(y, "nmarker")
  marker <- attr(y, "marker")
  
  # ----------------------------------------------------------------------------------------
  ## Set alpha/mu/sigma intercept starting value.
  # Noch mal mit Niki absprechen
  # Warum initialisiert man alpha mit machine.eps? Man muss dann bei MJM 
  # wahrscheinlich auch markerm2 initialisieren?
  # Computational besser für Ableitungen, so lassen
  eta_timegrid_alpha <- 0
  #eta_T_alpha <- 0
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
      #eta_T_sj <- drop(x$alpha$smooth.construct[[j]]$XT %*% b)
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- eta_grid_sj
      #x$alpha$smooth.construct[[j]]$state$fitted_T <- eta_T_sj
      eta_timegrid_alpha <- eta_timegrid_alpha + eta_grid_sj
      #eta_T_alpha <- eta_T_alpha + eta_T_sj
    }
  } 
  
  # Hier muss man sich noch überlegen, was man dann mit dem pcre-Term machen 
  # möchte: Wie wird der dann upgedated? Das Xgrid besteht aus col = (id, wfpc1,
  # wfpc2) und nrow = 700, und b = 2*50 - 2 Parameter
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
        # Unterschied smooth.construct, smoothCon Funktion mit Constraint
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
  iter <- 0
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
    # ---- LIKELIHOOD für longitudinal-Teil ist noch falsch berechnet!
    # Man verwendet hier eigentlich das eta_timegrid
    eta_T_long <- drop(
      t(rep(1, nmarker)) %x% diag(nsubj) %*% (eta$alpha*eta_T_mu))
    eta_T <- eta$lambda + eta$gamma + eta_T_long
    # sum_Lambda muss doch auch noch für die jeweilige Intervall-Länge gewichtet
    # werden, oder?
    sum_Lambda <- (survtime/2 * exp(eta$gamma)) %*%
      (diag(nsubj)%x%t(gq_weights))%*%
      exp(eta_timegrid)
    logLik <- drop(status %*% eta_T - sum_Lambda) +
      sum(dnorm(y[[1]][, "obs"], mean = eta$mu, sd = exp(eta$sigma),
                log = TRUE))
    # Eigentlich sollte hier doch auch über die Log-Posterior das Max gebildet
    # werden? Die Score und Hesse-Funktionen beziehen nämlich schon die Prioris
    # mit ein
    

    # eta1_long fehlt und muss vielleicht noch überarbeitet werden? JM.R(915)?
    iter <- iter + 1
    eta1_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eta1_alpha <- matrix(eta$alpha, nrow = nsubj, ncol = nmarker)
    eta1 <- cbind(eta1_surv, eta0_alpha)
    eps0 <- mean(abs((eta1 - eta0) / eta1), na.rm = TRUE)
    #if (iter %% 5 == 0) {
      cat("It ", iter,", LogLik ", logLik, "\n")
    #}
    eta0 <- eta1
    eta0_surv <- eta1_surv
    eta0_alpha <- eta1_alpha
  }

  # Log-Posterior ausrechnen und ausgeben
  # get_logPost Funktion aus JM als Vorlage
  # log Likelihood reicht aus über eta
  # in bamlss:::simsurv kann man sich Cox-Modell auch simulieren lassen
  ## return(list("parameters" = par, "fitted.values" = eta))
  
  max_y <- max(y[[1]][, 1])
  pred_data <- data.frame(seq(0, max_y, length.out = 100))
  colnames(pred_data) <- x$lambda$smooth.construct[[1]]$term
  k <- x$lambda$smooth.construct[[1]]$bs.dim
  b_it <- x$lambda$smooth.construct[[1]]$state$parameters[-k]
  pred_mat <- PredictMat(x$lambda$smooth.construct[[1]], pred_data)
  plot(pred_data[, 1], pred_mat%*%b_it)
  
  
  stop("Basst.")
}

update_mjm_lambda <- function(x, y, nu, eta, eta_timegrid, survtime, ...)
{
  ## grid matrix -> x$Xgrid
  ## design matrix -> x$X
  ## penalty matrices -> x$S
  ## optimizer.R -> bfit_iwls() updating.
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  #tau2 <- bamlss::get.state(x, "tau2")

  int_i <- survint_gq(pred = "lambda", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_vec = x$Xgrid, weights = attr(y, "gq_weights"),
                      survtime = survtime)
  

  # Status from MJM_transform
  # XT from sm_time_transform_mjm()
  x_score <- drop(attr(y, "status") %*% x$XT) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  ## Newton-Raphson.
  # Minus? z.B. in zeile 1211 g + nu * HS
  # bamlss::JM verwendet matrix_inv() Funktion definiert in BAMLSS.R
  # Ausgleich über nu?
  # -
  # Prior aus xhess verwenden
  # Dann doch wieder mit plus
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  return(x$state)
  
}

update_mjm_gamma <- function(x, y, nu, eta, eta_timegrid, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  take_last <- attr(y, "take_last")
  exp_eta_gamma <- exp(eta$gamma)
  
  int_i <- survint_gq(pred = "gamma", pre_fac = exp_eta_gamma, pre_vec = x$X,
                      omega = exp(eta_timegrid),
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  x_score <- drop(attr(y, "status") %*% x$X) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted.values <- drop(x$X %*% b)
  return(x$state)
  
}

update_mjm_alpha <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_mu, 
                             eta_T_mu, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  # int_i <- survint_gq(pred = "long", pre_fac = rep(exp(eta$gamma), nmarker), 
  #                     omega = rep(exp(eta_timegrid), nmarker),
  #                     int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
  #                     weights = attr(y, "gq_weights"),
  #                     survtime = survtime)
  int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
                      

  delta <- rep(attr(y, "status"), nmarker)
  # Verwende hier x$X und nicht x$XT?
  
  x_score <- drop(t(delta * x$XT) %*% eta_T_mu) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  # Hier wird in jm_bamlss auch x$ und nicht x$XT verwendet
  x$state$fitted.values <- drop(x$X %*% b)
  # Braucht man das fitted_T überhaupt? Sind nicht die fitted.values eh schon
  # die fitted_T - Werte, weil nämlich x$XT überhaupt nicht nötig war zu
  # erstellen, weil das eh schon in x$X enthalten war?
  # x$state$fitted_T <- drop(x$XT %*% b)

  return(x$state)
  
}

update_mjm_mu <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_alpha, 
                          survtime, ...) {

  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)

  delta <- rep(attr(y, "status"), nmarker)
  x_score <- drop(
    crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
      t(delta * x$XT) %*% eta$alpha) - colSums(int_i$score_int)
  x_H <- crossprod(x$X * (1 / exp(eta$sigma)^2), x$X) +
    matrix(colSums(int_i$hess_int), ncol = b_p)

  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta

  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$fitted_T <- drop(x$XT %*% b)

  return(x$state)
  
}

# Muss man hier nicht noch die Summe dafür anpassen, dass man nicht bis -1,
# sondern bis T_i integriert? Hier noch ein Argument einfügen, dass als Faktor
# noch T_i/2 auf jedes Element hinzumultipliziert?
survint_gq <- function(pred = c("lambda", "gamma", "long"), pre_fac,
                       pre_vec = NULL, omega, int_fac = NULL,
                       int_vec = NULL, weights, survtime) {
   
  if (sum(c(is.null(pre_vec), is.null(int_vec))) != 1) {
    stop("Either pre_vec or int_vec must be specified.")
  }
  if (pred != "long" & !is.null(int_fac)) {
    stop("Argument int_fac is only used for longitudinal predictors.")
  }
  # int_vec <- (t(rep(1, nmarker)) %x% diag(length(eta_timegrid))) %*%
  #   (eta_timegrid_alpha * x$Xgrid)
  # Integration as weighted sum of evaluated points
  n <- length(survtime)
  gq_mat <- diag(n)%x%t(weights)
  
  if (pred == "long") {
    nmarker <- nrow(int_vec)/(n*length(weights))
    if (nmarker %% 1 != 0) {
      stop("Dimensions of longitudinal design matrix do not match.")
    }
    pre_fac <- rep(pre_fac, nmarker)
    omega <- rep(omega, nmarker)
    gq_mat <- diag(n*nmarker)%x%t(weights)
    
    # Alternativ mit Matrixmultiplikation statt verlängertem Vektor
    # dim_mat <- t(rep(1, nmarker)) %x% diag(n*length(weights))
    # score_int <- pre_fac*gq_mat %*% (omega*dim_mat %*% (int_fac*int_vec))
    # hess_int <- pre_fac*gq_mat %*% 
    #   (omega*dim_mat %*% t(apply(int_fac*int_vec, 1, tcrossprod)))
  }
  switch(pred, 
    "lambda" =, "long" = {
      if (is.null(int_fac)) {
        int_fac <- rep(1, n*length(weights))
      }
      score_int <- pre_fac*gq_mat %*% (omega*int_fac*int_vec)
      hess_int <- pre_fac*gq_mat %*% (omega*int_fac^2*t(apply(int_vec, 1,
                                                              tcrossprod)))
    },
    "gamma" = {
      pre_fac <- c(pre_fac*gq_mat %*% omega)
      score_int <- pre_fac*pre_vec
      hess_int <- pre_fac*t(apply(pre_vec, 1, tcrossprod))
    })
  
  if (dim(score_int)[1] != dim(hess_int)[1]) {
    hess_int <- t(hess_int)
    if (dim(score_int)[1] != dim(hess_int)[1]) {
      stop("Problem with dimensions in gauss quadrature.")
    }
  }
  list(score_int = survtime/2 * score_int, hess_int = survtime/2 * hess_int)
  # hess_int: each row corresponds to one individual 
  # each row has ncol(vec)^2 elements -> is a matrix of derivatives
}

Surv2 <- bamlss:::Surv2

f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1, # marker does not make sense here
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)


# Externe zeitvariierende Kovariablen nicht möglich.
# in example:
# eta_lambda = 1.4*log((t + 10)/1000)
curve(1.4*log((x + 10)/1000), from = 0, to = 120)
b <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime")




if(FALSE) {
  x <- seq(0, 3, length = 100)
  y <- sin(x)

  gq <- statmod::gauss.quad(100, kind = "legendre")

  3/2 * sum(sin(3/2 * gq$nodes + 3/2) * gq$weights)

  b <- bamlss(y ~ s(x))

  foo <- function(x) {
    X <- cbind(PredictMat(b$x$mu$smooth.construct[["s(x)"]],
                          data.frame("x" = x)), 1)
    beta <- coef(b, FUN = mean, hyper = FALSE, model = "mu")
    return(drop(X %*% beta))
  }

  3/2 * sum(foo(3/2 * gq$nodes + 3/2) * gq$weights)
}

if(FALSE) {
  pref <- c(1, 1, 1)
  prev <- matrix(1:4, byrow = TRUE, ncol = 4, nrow = 3)
  om <- c(rep(1, 7), rep(2, 7), rep(3, 7))
  intf <- rep(1, 21)
  intv <- matrix(rep(1:4, each = 21), ncol = 4, nrow = 21)
  we <- rep(1, 7)#statmod::gauss.quad(7)$weights
  sur <- c(1, 1, 1)
  
  # lambda case
  lc <- survint_gq(pred = "lambda", pre_fac = pref, pre_vec = NULL, omega = om,
                   int_fac = NULL, int_vec = intv, weights = we, survtime = sur)
  
  # alpha case
  ac <- survint_gq(pred = "long", pre_fac = pref, pre_vec = NULL, omega = om,
                   int_fac = intf, int_vec = intv, weights = we, survtime = sur)
  
  # gamma case
  gc <- survint_gq(pred = "gamma", pre_fac = pref, pre_vec = prev, omega = om,
                   int_fac = NULL, int_vec = NULL, weights = we, survtime = sur)
}

if(FALSE) {
  simSurv2 <- function (n = 300) 
  {
    X <- matrix(NA, nrow = n, ncol = 3)
    X[, 1] <- runif(n, -1, 1)
    X[, 2] <- runif(n, -3, 3)
    X[, 3] <- runif(n, -1, 1)
    cens_fct <- function(time, mean_cens) {
      censor_time <- rexp(n = length(time), rate = 1/mean_cens)
      event <- (time <= censor_time)
      t_obs <- apply(cbind(time, censor_time), 1, min)
      return(cbind(t_obs, event))
    }
    lambda <- function(time, x) {
      exp(log(time) + sin(time * 2))
    }
    d <- rSurvTime2(lambda, X, cens_fct, mean_cens = 5)
    return(d)
  }
  set.seed(123)
  n <- 300
  simpledata <- simSurv2(n)
  simpledata$id <- factor(seq_len(n))
  simplef <- list(
    Surv2(time, event, obs = x1) ~ -1 + s(time, k = 10),
    gamma ~ 1,
    mu ~ s(id, x2, x3, bs = "pcre", xt = list("mfpc" = MFPCA))
  )
  b <- bamlss(simplef, family = mjm_bamlss, data = simpledata, timevar = "time")
  curve(log(x) + sin(x * 2), from = 0, to = max(simpledata$time))
  simplef_cox <- list(
    Surv2(time, event, obs = x1) ~ -1 + s(time, k = 10),
    gamma ~ 1)
  b_c <- bamlss(simplef_cox, family = cox_bamlss, data = simpledata)
  plot(b_c, ask = FALSE)
}
