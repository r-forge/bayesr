library("refund")
library("bamlss")
library("MFPCA")

source("simMultiJM.R")
source("eval_mfun.R")

if(FALSE) {
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
MJM_transform <- function(object, subdivisions = 10, timevar = NULL, ...)
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
  grid <- function(time) {
    time / 2 * gq$nodes + time / 2
  }
  y2 <- cbind(object$y[[1]][, "time"], object$model.frame[[idvar]])
  colnames(y2) <- c("time", idvar)

  take_l <- !duplicated(y2[, 1:2])
  take_last_l <- !duplicated(y2, fromLast = TRUE)
  nsubj <- length(unique(y2[, idvar]))

  marker <- FALSE
  if(!is.null(object$model.frame$marker)) {
    marker <- TRUE
    y2 <- cbind(y2, "marker" = as.factor(object$model.frame$marker))
  }

  take <- !duplicated(y2)
  take_last <- !duplicated(y2, fromLast = TRUE)

  y2 <- y2[take, , drop = FALSE]
  nobs <- nrow(y2)
  grid <- lapply(y2[, "time"], grid)

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
              grid, yname, timevar_mu, if(i == "lambda") take_last_l else take_last,
              N = nsubj * subdivisions)
          } else {
            object$x[[i]]$smooth.construct[[j]] <- sm_time_transform_mjm(object$x[[i]]$smooth.construct[[j]],
              object$model.frame[, unique(c(xterm, yname, by, if(i == "mu") timevar_mu else timevar, idvar)), drop = FALSE],
              grid, yname, if(i == "mu") timevar_mu else timevar, if(i == "lambda") take_last_l else take_last)
          }
        }
      }
    }
  }

  ## Now linear part.
  for(i in c("lambda", "mu", "alpha")) {
    if(!is.null(object$x[[i]]$smooth.construct$model.matrix)) {
      object$x[[i]]$smooth.construct$model.matrix <- param_time_transform_mjm(object$x[[i]]$smooth.construct$model.matrix,
        bamlss:::drop.terms.bamlss(object$x[[i]]$terms, sterms = FALSE, keep.response = FALSE), object$model.frame, grid, yname, 
        if(i != "mu") timevar else timevar_mu, if(i == "lambda") take_last_l else take_last)
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


opt_MJM <- function(x, y, eps = 0.0001, maxit = 400, nu = 0.1, ...)
{
  y <- y[[1]]

  eta <- bamlss:::get.eta(x, expand = FALSE)

  eps0 <- eps + 1
  eta0 <- do.call("cbind", eta)
  iter <- 0

  while((eps0 > eps) & (iter < maxit)) {
    ## (1) update lambda.
    for(j in names(x$lambda$smooth.construct)) {
      state <- update_mjm_lambda(x$lambda$smooth.construct[[j]], nu = nu, ...)
    }
  }
}

update_mjm_lambda <- function(x, nu, ...)
{
  ## grid matrix -> x$Xgrid
  ## design matrix -> x$X
  ## penalty matrices -> x$S
  ## optimizer.R -> bfit_iwls() updating.
  b <- bamlss::get.state(x, "b")
  tau2 <- bamlss::get.state(x, "tau2")

print(x$state)
print(x$term)
print(x$label)
print(tau2)
stop("!\n")

  ## Newton-Raphson.
  ## b <- b + nu * H %*% grad
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
