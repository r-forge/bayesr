
# Updating lambda predictor -----------------------------------------------


update_mjm_lambda <- function(x, y, nu, eta, eta_timegrid, survtime, ...)
{
  ## grid matrix -> x$Xgrid
  ## design matrix -> x$X
  ## penalty matrices -> x$S
  ## optimizer.R -> bfit_iwls() updating.
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  
  int_i <- survint_C(pred = "lambda", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_vec = x$Xgrid, weights = attr(y, "gq_weights"),
                      survtime = survtime)
  
  
  # Status from MJM_transform
  # XT from sm_time_transform_mjm()
  x_score <- drop(attr(y, "status") %*% x$XT) - int_i$score_int
  x_H <- matrix(int_i$hess_int, ncol = b_p)
  
  ## Newton-Raphson.
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  return(x$state)
  
}



# Updating gamm predictor -------------------------------------------------


update_mjm_gamma <- function(x, y, nu, eta, eta_timegrid, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  take_last <- attr(y, "take_last")
  exp_eta_gamma <- exp(eta$gamma)
  
  int_i <- survint_C(pred = "gamma", pre_fac = exp_eta_gamma, pre_vec = x$X,
                      omega = exp(eta_timegrid),
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  x_score <- drop(attr(y, "status") %*% x$X) - int_i$score_int
  x_H <- matrix(int_i$hess_int, ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted.values <- drop(x$X %*% b)
  return(x$state)
  
}


# Updating alpha predictor ------------------------------------------------


update_mjm_alpha <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_mu, 
                             eta_T_mu, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  
  
  delta <- rep(attr(y, "status"), nmarker)
  
  x_score <- drop(t(delta * x$XT) %*% eta_T_mu) - int_i$score_int
  x_H <- matrix(int_i$hess_int, ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  
  return(x$state)
  
}



# Updating mu predictor ---------------------------------------------------


update_mjm_mu <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_alpha, 
                          survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  
  delta <- rep(attr(y, "status"), nmarker)
  x_score <- drop(
    crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
      t(delta * x$XT) %*% eta$alpha) - int_i$score_int
  x_H <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
    matrix(int_i$hess_int, ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  S <- bamlss:::matrix_inv(x_H, index = NULL)
  b <- b + nu * S %*% x_score
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  x$state$fitted_T <- drop(x$XT %*% b)
  
  return(x$state)
  
}



# Updating sigma predictor ------------------------------------------------


update_mjm_sigma <- function(x, y, nu, eta, eta_timegrid, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  
  x_score <- crossprod(x$X, -1 + (y[[1]][, "obs"] - eta$mu)^2 / 
                         exp(eta$sigma)^2)
  x_H <- 2 * crossprod(x$X * drop((y[[1]][, "obs"] - eta$mu)/ exp(eta$sigma)^2),
                       x$X * drop(y[[1]][, "obs"] - eta$mu))

  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted.values <- drop(x$X %*% b)

  return(x$state)
  
}
