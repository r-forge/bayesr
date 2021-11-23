
# Survival Integral Function survint --------------------------------------

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
           hess_int <- pre_fac*gq_mat %*% sum_mat %*% 
             (int_fac^2*t(apply(int_vec, 1, tcrossprod)))
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




# Create simple univariate example for hazard function -------------------

n <- 3
n_gq <- 7
subdivisions <- 25
nmarker <- 1

# Survival times
s_times <- seq(0, 1, length.out = (n+1))[-1]

# Integration GQ
gq <- statmod::gauss.quad(n = n_gq)
create_grid <- function(time) {
  time / 2 * gq$nodes + time / 2
}
int_times <- unlist(lapply(s_times, create_grid))


# constant baseline hazard, constant association
lambda <- Vectorize(function (t) {
  1
})
gamma <- function (x) {
  1*x
}
alpha <- Vectorize(function (x) {
  1
})
mu <- function (t) {
  1 + 1*t
}

# This gives the function to integrate:
# Intercept: \int exp(3+x) dx -> exp(3+x)
# Time: \int x*exp(3+x) dx -> exp(x+3)*(x-1)
# => Hess: exp(x+3), exp(x+3)*(x-1), exp(x+3)*(x-1), exp(x+3)*(x^2-2x+2)

eta_timegrid_lambda <- lambda(rep(1, n*n_gq))
eta_gamma <- gamma(rep(1, n))
eta_timegrid_alpha <- alpha(rep(1, n*n_gq))
eta_timegrid_mu <- mu(int_times)
eta_timegrid_long <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu))
eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long

x_Xgrid <- cbind(1, int_times)

int_new <- survint_gq(pred = "long", 
                      pre_fac = exp(eta_gamma), # eta$gamma
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha, 
                      int_vec = x_Xgrid,
                      weights = gq$weights,
                      survtime = s_times)

# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(exp(3+t)-exp(3), exp(3+t)*(t-1) + exp(3)),
       hess_int = c(exp(3+t)-exp(3), exp(3+t)*(t-1) + exp(3),
                    exp(3+t)*(t-1) + exp(3),
                    exp(3+t)*(t^2-2*t+2) - 2*exp(3)))
}
int_func(s_times[1])



# Implementation in bamlss

# Integration grid
int_times_bamlss <- rep(1, n*subdivisions)
for (i in seq_len(n)) {
  int_times_bamlss[((i-1)*subdivisions + 1):(i*subdivisions)] <- 
    seq(0, s_times[i], length.out = subdivisions)
}

eta_timegrid_lambda_bamlss <- lambda(rep(1, n*subdivisions))
eta_timegrid_alpha_bamlss <- alpha(rep(1, n*subdivisions))
eta_timegrid_mu_bamlss <- mu(int_times_bamlss)
eta_timegrid_bamlss <- eta_timegrid_lambda_bamlss +
  eta_timegrid_alpha_bamlss * eta_timegrid_mu_bamlss

X <- cbind(1, int_times_bamlss)
dimnames(X) <- list(NULL, c("Intercept", "Time"))
width <- int_times_bamlss[(seq_len(n) - 1)*subdivisions + 2]
ind <- matrix(rep(as.integer(c(1, -1), n)), ncol = 2, nrow = n, byrow = TRUE)

# For some reason, this does not work
# int_bamlss <- bamlss:::survint(X = X, 
#                                eta = exp(eta_timegrid_bamlss) *
#                                  eta_timegrid_alpha_bamlss, 
#                                width = width,
#                                gamma = exp(eta_gamma),
#                                eta2 = exp(eta_timegrid_bamlss) * 
#                                  eta_timegrid_alpha_bamlss^2,
#                                index = ind,
#                                dX = NULL, 
#                                deta = NULL,
#                                deta2 = NULL)


# Now Multiple Markers ----------------------------------------------------

nmarker <- 2

# Have the same function for the second longitudinal marker and the association
# This gives the function to integrate for marker-specific effects:
# Intercept: \int exp(4+2x) dx -> 0.5*exp(4+2x)
# Time: \int x*exp(4+2x) dx -> 0.25*exp(2x+4)*(2x-1)
# => Hess: 0.5*exp(4+2x), 0.25*exp(2x+4)*(2x-1), 0.25*exp(2x+4)*(2x-1),
#          0.25*exp(2x+4)*(2x^2-2x+1)

eta_timegrid_alpha <- alpha(rep(rep(1, n*n_gq), nmarker))
eta_timegrid_mu <- mu(rep(int_times, nmarker))
eta_timegrid_long <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu))
eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long

x_Xgrid <- rbind(cbind(1, int_times), cbind(0, rep(0, length(int_times))))

int_new <- survint_gq(pred = "long", 
                      pre_fac = exp(eta_gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha, 
                      int_vec = x_Xgrid,
                      weights = gq$weights,
                      survtime = s_times)

# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(0.5*exp(4+2*t)-0.5*exp(4), 
                     0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4)),
       hess_int = c(0.5*exp(4+2*t)-0.5*exp(4), 
                    0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4),
                    0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4),
                    0.25*exp(2*t+4)*(2*t^2-2*t+1) - 0.25*exp(4)))
}
int_func(s_times[1])


# This gives the function to integrate for marker-common effects:
# Intercept: \int 2*exp(4+2x) dx -> exp(4+2x)
# Time: \int 2x*exp(4+2x) dx -> 0.5*exp(2x+4)*(2x-1)
# => Hess: exp(4+2x), 0.5*exp(2x+4)*(2x-1), 0.5*exp(2x+4)*(2x-1),
#          0.5*exp(2x+4)*(2x^2-2x+1)

x_Xgrid_re <- cbind(1, rep(int_times, nmarker))
int_new_re <- survint_gq(pred = "long", 
                      pre_fac = exp(eta_gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha, 
                      int_vec = x_Xgrid_re,
                      weights = gq$weights,
                      survtime = s_times)


# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(exp(4+2*t)-exp(4), 0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4)),
       hess_int = c(exp(4+2*t)-exp(4), 0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4),
                    0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4),
                    0.5*exp(2*t+4)*(2*t^2-2*t+1) - 0.5*exp(4)))
}
int_func(s_times[1])
