
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
sm_time_transform_mjm <- function(x, data, grid, yname, timevar, take)
{

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
  }
  if(x$by != "NA" & x$by != yname)
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))

  x$Xgrid <- PredictMat(x, X)
  
  # XT necessary for calculation of score and hessian
  # Matrix of evaluations at the vector of survival times
  x$XT <- PredictMat(x, data)
  
  x
}

## PCRE transformer.
sm_time_transform_mjm_pcre <- function(x, data, grid, yname, timevar, take, N)
{
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
  }
  if(x$by != "NA" & x$by != yname)
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))

  class(x) <- "pcre2.random.effect"

  x$N <- N
  x$term <- c(x$term, timevar)
  x$timevar <- timevar
  x$Xgrid <- PredictMat(x, X)

  x
}

Predict.matrix.pcre2.random.effect <- function(object, data)
{
  if(is.null(object$xt$mfpc))
    stop("need mfpa object!")
  X <- eval_mfpc(object$xt$mfpc, data[[object$timevar]][1:object$N])
  if(ncol(X) != (length(object$term) - 2))
    stop("check M argument in MFPCA()!")
  X <- cbind(data[[object$term[1]]], X)
  colnames(X) <- object$term[-length(object$term)]
  return(X)
}

## Linear design transformer.
param_time_transform_mjm <- function(x, formula, data, grid, yname, timevar, take)
{
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- Xn <- NULL
  for(j in names(data)) {
    if((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", j, fixed = TRUE)) & (j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
      Xn <- c(Xn, j)
    }
  }
  if(!is.null(X))
    colnames(X) <- Xn
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }

  x$Xgrid <- model.matrix(formula, data = X)

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
  

  marker <- FALSE
  if(!is.null(object$model.frame$marker)) {
    y2_l <- cbind(y2_l, "marker" = as.factor(object$model.frame$marker))
    marker <- length(levels(as.factor(object$model.frame$marker)))
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
  attr(object$y, "marker") <- marker
  

  yname <- all.names(object$x$lambda$formula[2])[2]
  timevar_mu <- timevar
  if(is.null(timevar_mu))
    stop("the time variable is not specified, needed for mu!")
  timevar <- yname

  ## The basic setup.
  if(is.null(attr(object$x, "bamlss.engine.setup")))
    object$x <- bamlss.engine.setup(object$x, ...)

  ## Remove intercept from lambda.
  if(!is.null(object$x$lambda$smooth.construct$model.matrix)) {
    cn <- colnames(object$x$lambda$smooth.construct$model.matrix$X)
    if("(Intercept)" %in% cn)
      object$x$lambda$smooth.construct$model.matrix$X <- object$x$lambda$smooth.construct$model.matrix$X[, cn != "(Intercept)", drop = FALSE]
    if(ncol(x$lambda$smooth.construct$model.matrix$X) < 1) {
      object$x$lambda$smooth.construct$model.matrix <- NULL
      object$x$lambda$terms <- drop.terms.bamlss(x$lambda$terms, pterms = FALSE, keep.intercept = FALSE)
    } else {
      object$x$lambda$smooth.construct$model.matrix$term <- gsub("(Intercept)+", "",
        object$x$lambda$smooth.construct$model.matrix$term, fixed = TRUE)
      object$x$lambda$smooth.construct$model.matrix$state$parameters <- object$x$lambda$smooth.construct$model.matrix$state$parameters[-1]
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
          if(inherits(object$x[[i]]$smooth.construct[[j]], "pcre.random.effect")) {
            object$x[[i]]$smooth.construct[[j]] <- sm_time_transform_mjm_pcre(object$x[[i]]$smooth.construct[[j]],
              object$model.frame[, unique(c(xterm, yname, by, timevar_mu, idvar)), drop = FALSE],
              grid_l, yname, timevar_mu, take_last_l, N = nsubj * subdivisions)
          } else {
            object$x[[i]]$smooth.construct[[j]] <- sm_time_transform_mjm(object$x[[i]]$smooth.construct[[j]],
              object$model.frame[, unique(c(xterm, yname, by, if(i == "mu") timevar_mu else timevar, idvar)), drop = FALSE],
              if(i == "lambda") grid else grid_l, yname, 
              if(i == "mu") timevar_mu else timevar,
              if(i == "lambda") take_last else take_last_l)
          }
        }
      }
    }
  }

  ## Now linear part.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct$model.matrix)) {
      object$x[[i]]$smooth.construct$model.matrix <- param_time_transform_mjm(object$x[[i]]$smooth.construct$model.matrix,
        bamlss:::drop.terms.bamlss(object$x[[i]]$terms, sterms = FALSE, keep.response = FALSE), object$model.frame,
        if(i == "lambda") grid else grid_l, yname, 
        if(i != "mu") timevar else timevar_mu, if(i == "lambda") take_last else take_last_l)
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


opt_MJM <- function(x, y, start = NULL, eps = 0.0001, maxit = 400, nu = 0.1, ...)
{
  
  if(!is.null(start))
    x <- bamlss:::set.starting.values(x, start)
  
  # ----------------------------------------------------------------------------------------
  ## Set alpha/mu/sigma intercept starting value.
  # Noch mal mit Niki absprechen
  # Warum initialisiert man alpha mit machine.eps? Man muss dann bei MJM 
  # wahrscheinlich auch markerm2 initialisieren?
  # Computational besser für Ableitungen, so lassen
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$alpha$smooth.construct[[j]]$state$parameters[1] <- .Machine$double.eps
        x$alpha$smooth.construct[[j]]$state$fitted.values <- 
          x$alpha$smooth.construct[[j]]$X %*% 
          x$alpha$smooth.construct[[j]]$state$parameters
      } 
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      eta_sj <- drop(x$alpha$smooth.construct[[j]]$Xgrid %*% b)
      x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- eta_sj
      eta_timegrid_alpha <- eta_timegrid_alpha + eta_sj
    }
  } 
  
  # Man braucht hier noch die Info über die Marker. Dann kann man jeden Marker 
  # Intercept mit seinem eigenen Mittelwert initialisieren
  #Initialisiserung kann man auch weglassen - Input sollte standardisiert sein
  
  # Hier muss man sich noch überlegen, was man dann mit dem pcre-Term machen 
  # möchte: Wie wird der dann upgedated? Das Xgrid besteht aus col = (id, wfpc1,
  # wfpc2) und nrow = 700, und b = 2*50 - 2 Parameter
  eta_timegrid_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      if (j == "model.matrix" & is.null(start)) {
        x$mu$smooth.construct[[j]]$state$parameters[1] <- 
          mean(y[[1]][, "obs"], na.rm = TRUE)
        x$mu$smooth.construct[[j]]$state$fitted.values <- 
          x$mu$smooth.construct[[j]]$X %*% 
          x$mu$smooth.construct[[j]]$state$parameters
      }
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      if(inherits(x$mu$smooth.construct[[j]], "pcre2.random.effect")){
        # Unterschied smooth.construct, smoothCon Funktion mit Constraint
        eta_sj <- rep(0, nrow(x$mu$smooth.construct[[j]]$Xgrid))
      } else {
        eta_sj <- drop(x$mu$smooth.construct[[j]]$Xgrid %*% b)
      }
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- eta_sj
      eta_timegrid_mu <- eta_timegrid_mu + eta_sj
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

  # Auch hier sollte man für jeden Marker einen eigenen Startwert wählen
  if(!is.null(x$sigma$smooth.construct$model.matrix)) {
    x$sigma$smooth.construct$model.matrix$state$parameters[1] <- 
      log(sd(y[[1]][, "obs"], na.rm = TRUE))
    x$sigma$smooth.construct$model.matrix$state$fitted.values <-
      x$sigma$smooth.construct$model.matrix$X %*% 
      x$sigma$smooth.construct$model.matrix$state$parameters
  }

  
  # Das eta sollte doch eigentlich unterschiedliche Anzahl an Werten haben?
  # Also für die survival-Parts nur die Anzahl der Beobachtungen und für den 
  # longitudinalen Teil die Gesamtzahl?
  # mjm_transform Zeile 200
  # Vor dem Entfernen des Intercepts -> Überprüfen. Sollte passen
  eta <- bamlss:::get.eta(x, expand = FALSE)

  marker <- attr(y, "marker")
  eta_timegrid_long <- eta_timegrid_alpha*eta_timegrid_mu
  if (marker) {
    eta_timegrid_long <- drop(
      t(rep(1, marker)) %x% diag(length(eta_timegrid_lambda)) %*%
      eta_timegrid_long)
  }
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
  
  eps0 <- eps + 1
  eta0 <- do.call("cbind", eta)
  iter <- 0
  
  # Für logLik
  nsubj <- attr(y, "nsubj")
  gq_weights <- attr(y, "gq_weights")
  take_last <- attr(y, "take_last")
  take_last_l <- attr(y, "take_last_l")
  status <- attr(y, "status")
  
  
  while((eps0 > eps) & (iter < maxit)) {
    ## (1) update lambda.
    for(j in names(x$lambda$smooth.construct)) {
      state <- update_mjm_lambda(x$lambda$smooth.construct[[j]], y = y, nu = nu, 
                                 eta = eta, eta_timegrid = eta_timegrid, ...)
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
                                 eta = eta, eta_timegrid = eta_timegrid, ...)
        eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[j]]$state) +
          fitted(state)
        x$gamma$smooth.construct[[j]]$state <- state
      }
    }
    
    
    # Likelihood calculation
    eta_T_long <- eta$alpha[take_last_l] * eta$mu[take_last_l]
    if (marker) {
      eta_T_long <- drop(
        t(rep(1, marker)) %x% diag(nsubj) %*% eta_T_long)
    }
    eta_T <- eta$lambda[take_last] + eta$gamma[take_last] + eta_T_long
    sum_Lambda <- exp(eta$gamma[take_last])%*%(diag(nsubj)%x%t(gq_weights))%*%
      exp(eta_timegrid)
    logLik <- status %*% eta_T - sum_Lambda # + longitudinal part
    # logLik <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    #   exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
    
    iter <- iter + 1
    eta1 <- do.call("cbind", eta)
    eps0 <- mean(abs((eta1 - eta0) / eta1), na.rm = TRUE)
    cat("It ", iter,", LogLik ", logLik, "\n")
    eta0 <- eta1
  }

  # Log-Posterior ausrechnen und ausgeben
  # get_logPost Funktion aus JM als Vorlage
  # log Likelihood reicht aus über eta
  # in bamlss:::simsurv kann man sich Cox-Modell auch simulieren lassen
  ## return(list("parameters" = par, "fitted.values" = eta))
  
  stop("Basst.")
}

update_mjm_lambda <- function(x, y, nu, eta, eta_timegrid, ...)
{
  ## grid matrix -> x$Xgrid
  ## design matrix -> x$X
  ## penalty matrices -> x$S
  ## optimizer.R -> bfit_iwls() updating.
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  tau2 <- bamlss::get.state(x, "tau2")
  
  exp_eta_gamma <- exp(eta$gamma[attr(y, "take_last")])
  
  int_i <- survint_gq(pre_fac = exp_eta_gamma, omega = exp(eta_timegrid),
                      int_vec = x$Xgrid, weights = attr(y, "gq_weights"))
  # For other predictors
  # * to be created
  # -----
  # mu:
  # int_i <- survint_gq(pre_fac = exp_eta_gamma, omega = exp(eta_timegrid),
  #                     int_fac = eta_alpha*, int_vec = x$Xgrid,
  #                     weights = attr(y, "gq_weights"))
  # -----
  # alpha:
  # int_i <- survint_gq(pre_fac = exp_eta_gamma, omega = exp(eta_timegrid),
  #                     int_fac = eta_mu*, int_vec = x$Xgrid, 
  #                     weights = attr(y, "gq_weights"))

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

update_mjm_gamma <- function(x, y, nu, eta, eta_timegrid, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  take_last <- attr(y, "take_last")
  exp_eta_gamma <- exp(eta$gamma[take_last])
  x_gamma <- x$X[take_last, , drop = FALSE]
  
  int_i <- survint_gq(pre_fac = exp_eta_gamma, pre_vec = x_gamma,
                      omega = exp(eta_timegrid),
                      weights = attr(y, "gq_weights"))
  x_score <- drop(attr(y, "status") %*% x_gamma) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted.values <- drop(x$X %*% b)
  return(x$state)
  
}

survint_gq <- function(pre_fac, pre_vec = NULL, omega, int_fac = NULL,
                       int_vec = NULL, weights) {
  
  if (sum(c(is.null(pre_vec), is.null(int_vec))) != 1) {
    stop("Either pre_vec or int_vec must be specified.")
  }
  
  # Integration as weighted sum of evaluated points
  n <- length(pre_fac)
  gq_mat <- diag(n)%x%t(weights)
  
  if (is.null(pre_vec)) {
    # Not a gamma predictor
    if (is.null(int_fac)) {
      int_fac <- rep(1, n*length(weights))
    }
    score_int <- pre_fac*gq_mat %*% (omega*int_fac*int_vec)
    hess_int <- pre_fac*gq_mat %*% (omega*int_fac^2*t(apply(int_vec, 1, tcrossprod)))
  } else {
    # Gamma predictor
    if (!is.null(int_fac)) {
      stop("Argument int_fac is not used for predictor gamma.")
    }
    pre_fac <- c(pre_fac*gq_mat %*% omega)
    score_int <- pre_fac*pre_vec
    hess_int <- pre_fac*t(apply(pre_vec, 1, tcrossprod))
  }
  if (dim(score_int)[1] != dim(hess_int)[1]) {
    hess_int <- t(hess_int)
    if (dim(score_int)[1] != dim(hess_int)[1]) {
      stop("Problem with dimensions in gauss quadrature.")
    }
  }
  list(score_int = score_int, hess_int = hess_int)
  # hess_int: each row corresponds to one individual 
  # each row has ncol(vec)^2 elements -> is a matrix of derivatives
}

Surv2 <- bamlss:::Surv2

f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA)),
  sigma ~ 1,
  alpha ~ -1 + marker + s(survtime, by = marker)
)


# Externe zeitvariierende Kovariablen nicht möglich.

b <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime")




if(FALSE) {
  x <- seq(0, 3, length = 100)
  y <- sin(x)

  gq <- statmod::gauss.quad(100, kind = "legendre")

  3/2 * sum(sin(3/2 * gq$nodes + 3/2) * gq$weights)

  b <- bamlss(y ~ s(x))

  foo <- function(x) {
    X <- cbind(PredictMat(b$x$mu$smooth.construct[["s(x)"]], data.frame("x" = x)), 1)
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
  
  # lambda case
  lc <- survint_gq(pre_fac = pref, pre_vec = NULL, omega = om,
                   int_fac = NULL, int_vec = intv, weights = we)
  
  # alpha case
  ac <- survint_gq(pre_fac = pref, pre_vec = NULL, omega = om,
                   int_fac = intf, int_vec = intv, weights = we)
  
  # gamma case
  gc <- survint_gq(pre_fac = pref, pre_vec = prev, omega = om,
                   int_fac = NULL, int_vec = NULL, weights = we)
}

if(FALSE) {
  set.seed(123)
  simpledata <- simSurv()
  simpledata$id <- factor(1:300)
  simplef <- list(
    Surv2(time, event, obs = x1) ~ -1 + s(time),
    gamma ~ 1,
    mu ~ s(id, x2, x3, bs = "pcre", xt = list("mfpc" = MFPCA))
  )
  # Kann man dann außerhalb der Schleife machen
  b <- bamlss(simplef, family = mjm_bamlss, data = simpledata, timevar = "time")
}
