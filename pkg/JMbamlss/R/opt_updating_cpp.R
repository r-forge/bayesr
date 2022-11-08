

# Updating mu predictor ---------------------------------------------------


update_mjm_mu_old <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_alpha, 
                          survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  if (any(class(x) == "unc_pcre.random.effect")) {
    int_i <- survint_C(pred = "fpc_re", pre_fac = exp(eta$gamma),
                       omega = exp(eta_timegrid),
                       int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                       weights = attr(y, "gq_weights"),
                       survtime = survtime)
    x_H <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
      diag(int_i$hess_int)
  } else {
    int_i <- survint_C(pred = "long", pre_fac = exp(eta$gamma),
                       omega = exp(eta_timegrid),
                       int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                       weights = attr(y, "gq_weights"),
                       survtime = survtime)
    x_H <- eigenMapMatMult(t(x$X * (1 / exp(eta$sigma)^2)), x$X) +
      matrix(int_i$hess_int, ncol = b_p)
  }
  
  
  delta <- rep(attr(y, "status"), nmarker)
  x_score <- drop(
    crossprod(x$X, (y[[1]][, "obs"] - eta$mu) / exp(eta$sigma)^2)  + 
      crossprod(delta * x$XT, eta$alpha)) - int_i$score_int

  
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
