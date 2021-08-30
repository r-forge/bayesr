library("refund")

source("simMultiJM.R")

if(FALSE) {
  d <- simMultiJM()
}

jm_bamlss <- function(...)
{
  links = c(
    lambda = "log",
    gamma = "log",
    mu = "identity",
    sigma = "log",
    alpha = "identity"
  )
  
  rval <- list(
    "family" = "jm",
    "names" = c("lambda", "gamma", "mu", "sigma", "alpha"),
    "links" = links,
    "transform" = JM_transform,
    "optimizer" = NULL,
    "sampler" = NULL,
    "predict" = NULL
  )
  
  class(rval) <- "family.bamlss"
  rval
}

JM_transform <- function(object, subdivisions = 25, ...)
{
  stopifnot(requireNamespace("statmod"))

  gq <- statmod::gauss.quad(subdivisions, ...)
  gq$weights <- gq$weights * exp(gq$nodes^2)

  ## Get idvar.
  class_mu <- sapply(object$x$mu$smooth.construct, class)
  class_mu <- sapply(class_mu, function(x) x[1])
  if(!any(class_mu == "pcre.random.effect"))
    stop("need pcre smooth!")
  j <- which(class_mu == "pcre.random.effect")
  smj <- object$x$mu$smooth.construct[[j]]
  idvar <- smj$term[1]

  ## Setup integration.
  grid <- function(upper, length){
    seq(from = 0, to = upper, length = length)
  }
  y2 <- cbind(object$y[[1]][, "time"], object$model.frame[[idvar]])
  colnames(y2) <- c("time", idvar)
  take <- !duplicated(y2)
  take_last <- !duplicated(y2, fromLast = TRUE)
  y2 <- y2[take, , drop = FALSE]
  nobs <- nrow(y2)
  grid <- lapply(y2[, "time"], grid, length = subdivisions)

  yname <- all.names(object$x$lambda$formula[2])[2]

  ## Transform lambda.
  ## FIXME: remove intercept!
  ## Compute new design matrices for integration.
  if(!is.null(object$x$lambda$smooth.construct)) {
    for(j in names(object$x$lambda$smooth.construct)) {
      if(j != "model.matrix") {
        xterm <- x$lambda$smooth.construct[[j]]$term
        by <- if(x$lambda$smooth.construct[[j]]$by != "NA") {
          x$lambda$smooth.construct[[j]]$by
        } else NULL
        object$x$lambda$smooth.construct[[j]] <- sm_time_transform(object$x$lambda$smooth.construct[[j]],
          data[, unique(c(xterm, yname, by, timevar, idvar)), drop = FALSE],
          grid, yname, timevar, take)
      }
    }
  }
print(grid)

stop()
}

Surv2 <- bamlss:::Surv2

f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) + s(id, wfpc.1, wfpc.2, bs = "pcre"),
  sigma ~ 1,
  alpha ~ 1
)

b <- bamlss(f, family = jm_bamlss, data = d)




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
