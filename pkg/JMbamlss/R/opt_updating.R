
# Updating lambda predictor -----------------------------------------------


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



# Updating gamm predictor -------------------------------------------------


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


# Updating alpha predictor ------------------------------------------------


update_mjm_alpha <- function(x, y, nu, eta, eta_timegrid, eta_timegrid_mu, 
                             eta_T_mu, survtime, ...) {
  
  b <- bamlss::get.state(x, "b")
  b_p <- length(b)
  nmarker <- attr(y, "nmarker")
  
  int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_mu, int_vec = x$Xgrid,
                      weights = attr(y, "gq_weights"),
                      survtime = survtime)
  
  
  delta <- rep(attr(y, "status"), nmarker)
  
  x_score <- drop(t(delta * x$XT) %*% eta_T_mu) - colSums(int_i$score_int)
  x_H <- matrix(colSums(int_i$hess_int), ncol = b_p)
  
  x_score <- x_score + x$grad(score = NULL, x$state$parameters, full = FALSE)
  x_H <- x_H + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  delta <- solve(x_H, x_score)
  b <- b + nu * delta
  
  x$state$parameters[seq_len(b_p)] <- b
  x$state$fitted_timegrid <- drop(x$Xgrid %*% b)
  x$state$fitted.values <- drop(x$X %*% b)
  # Braucht man das fitted_T überhaupt? Sind nicht die fitted.values eh schon
  # die fitted_T - Werte, weil nämlich x$XT überhaupt nicht nötig war zu
  # erstellen, weil das eh schon in x$X enthalten war?
  # x$state$fitted_T <- drop(x$XT %*% b)
  
  return(x$state)
  
}



# Updating mu predictor ---------------------------------------------------


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


# Survival integral for GQ ------------------------------------------------

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
  
  switch(pred, 
         "lambda" = {
           score_int <- pre_fac*gq_mat %*% (omega*int_vec)
           hess_int <- pre_fac*gq_mat %*% (omega*t(apply(int_vec, 1, 
                                                         tcrossprod)))
         },
         "gamma" = {
           pre_fac <- c(pre_fac*gq_mat %*% omega)
           score_int <- pre_fac*pre_vec
           hess_int <- pre_fac*t(apply(pre_vec, 1, tcrossprod))
         },
         "long" = {
           nmarker <- nrow(int_vec)/(n*length(weights))
           if (nmarker %% 1 != 0) {
             stop("Dimensions of longitudinal design matrix do not match.")
           }
           
           # Alternativ mit Matrixmultiplikation statt verlängertem Vektor
           # dim_mat <- t(rep(1, nmarker)) %x% diag(n*length(weights))
           # score_int <- pre_fac*gq_mat %*% (omega*dim_mat %*% (int_fac*int_vec))
           # hess_int <- pre_fac*gq_mat %*% 
           #   (omega*dim_mat %*% t(apply(int_fac*int_vec, 1, tcrossprod)))
           sum_mat <- matrix(rep(diag(omega), nmarker), nrow = length(omega))
           score_int <- pre_fac*gq_mat %*% sum_mat %*% (int_fac*int_vec)
           hess_int <- if (dim(int_vec)[2] == 1) {
             pre_fac*gq_mat %*% sum_mat %*% 
               (int_fac^2*do.call(rbind, apply(int_vec, 1, tcrossprod,
                                              simplify = FALSE)))
           } else {
             pre_fac*gq_mat %*% sum_mat %*% 
               (int_fac^2*t(apply(int_vec, 1, tcrossprod)))
           }

         })
  
  if (dim(score_int)[1] != dim(hess_int)[1]) {
    hess_int <- t(hess_int)
    if (dim(score_int)[1] != dim(hess_int)[1]) {
      stop("Problem with dimensions in gauss quadrature.")
    }
  }
  # Multiply with borders of integration for Legendre
  list(score_int = survtime/2 * score_int, hess_int = survtime/2 * hess_int)
  # hess_int: each row corresponds to one individual 
  # each row has ncol(vec)^2 elements -> is a matrix of derivatives
}

