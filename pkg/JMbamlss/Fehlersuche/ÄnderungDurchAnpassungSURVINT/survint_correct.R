
# Compare the different integral functions --------------------------------
source("R/opt_updating.R")

# Check the survint function ----------------------------------------------

# Setup for one longitudinal marker
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
# Correct integrals would be
int_func_one <- function (t) {
  list(score_int = c(exp(3+t)-exp(3), exp(3+t)*(t-1) + exp(3)),
       hess_int = c(exp(3+t)-exp(3), exp(3+t)*(t-1) + exp(3),
                    exp(3+t)*(t-1) + exp(3),
                    exp(3+t)*(t^2-2*t+2) - 2*exp(3)))
}

# Set up arguments for survint
eta_timegrid_lambda <- lambda(int_times)
eta_gamma <- gamma(rep(1, n))
eta_timegrid_alpha <- alpha(rep(1, n*n_gq))
eta_timegrid_mu <- mu(int_times)
eta_timegrid_long <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu))
eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
x_Xgrid <- cbind(1, int_times)

int_new <- survint_gq(pred = "long",
                      pre_fac = exp(eta_gamma),
                      omega = exp(eta_timegrid),
                      int_fac = eta_timegrid_alpha,
                      int_vec = x_Xgrid,
                      weights = gq$weights,
                      survtime = s_times)
int_new
# Compare to correct integrals
int_func_one(s_times[1])


# Check the survint function for non-unit interval ------------------------

# Survival times
s_times_120 <- seq(0, 120, length.out = (n+1))[-1]
int_times_120 <- unlist(lapply(s_times_120, create_grid))
eta_timegrid_mu_120 <- mu(int_times_120)
eta_timegrid_long_120 <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu_120))
eta_timegrid_120 <- eta_timegrid_lambda + eta_timegrid_long_120
x_Xgrid_120 <- cbind(1, int_times_120)

int_func_one(s_times_120[1])
survint_gq(pred = "long",
           pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid_120),
           int_fac = eta_timegrid_alpha,
                      int_vec = x_Xgrid_120,
                      weights = gq$weights,
                      survtime = s_times_120)

# ------LARGE DIFFERENCES--------



# Numerical reasons for large difference? ---------------------------------


# Check example with lower integral value to see if it depends on numerical
# reasons
# New example with different values for the predictors (decreasing hazard)

lambda_dec <- Vectorize(function (t) {
  -1
})
gamma_dec <- function (x) {
  -1*x
}
alpha_dec <- Vectorize(function (x) {
  1
})
mu_dec <- function (t) {
  -1 - 0.01*t
}

# This gives the function to integrate:
# Intercept: \int exp(-3-0.01x) dx -> -100exp(-3-0.01x)
# Time: \int x*exp(-3-x) dx -> -100exp(-0.01x-3)*(x+100)
# => Hess: -100exp(-0.01x-3), -100exp(-0.01x-3)*(x+100), 
#         -100exp(-0.01x-3)*(x+100), -100exp(-0.01x-3)*(x^2+200x+20000)
# Correct integrals would be
int_func_one_dec <- function (t) {
  list(score_int = c(-100*exp(-3-0.01*t)+100*exp(-3), 
                     -100*exp(-3-0.01*t)*(t+100) + 10000*exp(-3)),
       hess_int = c(-100*exp(-3-0.01*t)+100*exp(-3), 
                    -100*exp(-3-0.01*t)*(t+100) + 10000*exp(-3),
                    -100*exp(-3-0.01*t)*(t+100) + 10000*exp(-3),
                    -100*exp(-3-0.01*t)*(t^2+200*t+20000) + 2000000*exp(-3)))
}


# Set up arguments for survint
eta_timegrid_lambda_dec <- lambda_dec(int_times)
eta_gamma_dec <- gamma_dec(rep(1, n))
eta_timegrid_alpha_dec <- alpha_dec(rep(1, n*n_gq))
eta_timegrid_mu_dec <- mu_dec(int_times)
eta_timegrid_long_dec <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha_dec*eta_timegrid_mu_dec))
eta_timegrid_dec <- eta_timegrid_lambda_dec + eta_timegrid_long_dec

survint_gq(pred = "long",
           pre_fac = exp(eta_gamma_dec),
           omega = exp(eta_timegrid_dec),
           int_fac = eta_timegrid_alpha_dec,
           int_vec = x_Xgrid,
           weights = gq$weights,
           survtime = s_times)

int_func_one_dec(s_times[1])

# Now with non-unit interval
eta_timegrid_mu_dec_120 <- mu_dec(int_times_120)
eta_timegrid_long_dec_120 <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha_dec*eta_timegrid_mu_dec_120))
eta_timegrid_dec_120 <- eta_timegrid_lambda_dec + eta_timegrid_long_dec_120

int_func_one_dec(s_times_120[1])
survint_gq(pred = "long",
           pre_fac = exp(eta_gamma_dec),
           omega = exp(eta_timegrid_dec_120),
           int_fac = eta_timegrid_alpha_dec,
           int_vec = x_Xgrid_120,
           weights = gq$weights,
           survtime = s_times_120)

# ------NO DIFFERENCES--------
all.equal(survint_gq(pred = "long",
                     pre_fac = exp(eta_gamma_dec),
                     omega = exp(eta_timegrid_dec_120),
                     int_fac = eta_timegrid_alpha_dec,
                     int_vec = x_Xgrid_120,
                     weights = gq$weights,
                     survtime = s_times_120)$hess_int[1,],
          int_func_one_dec(s_times_120[1])$hess_int)


# Integration using bamlss Implementation ---------------------------------

# FOR SOME REASON THIS DOES NOT WORK - R CRASHES

# # Integration grid
# int_times_bamlss <- rep(1, n*subdivisions)
# for (i in seq_len(n)) {
#   int_times_bamlss[((i-1)*subdivisions + 1):(i*subdivisions)] <-
#     seq(0, s_times[i], length.out = subdivisions)
# }
# eta_timegrid_lambda_bamlss <- lambda(rep(1, n*subdivisions))
# eta_timegrid_alpha_bamlss <- alpha(rep(1, n*subdivisions))
# eta_timegrid_mu_bamlss <- mu(int_times_bamlss)
# eta_timegrid_bamlss <- eta_timegrid_lambda_bamlss +
#   eta_timegrid_alpha_bamlss * eta_timegrid_mu_bamlss
# X <- cbind(1, int_times_bamlss)
# dimnames(X) <- list(NULL, c("Intercept", "Time"))
# width <- int_times_bamlss[(seq_len(n) - 1)*subdivisions + 2]
# ind <- matrix(rep(as.integer(c(1, -1), n)), ncol = 2, nrow = n, byrow = TRUE)
# 
# int_old <- bamlss:::survint(X = X,
#                             eta = exp(eta_timegrid_bamlss) *
#                               eta_timegrid_alpha_bamlss,
#                             width = width,
#                             gamma = exp(eta_gamma),
#                             eta2 = exp(eta_timegrid_bamlss) *
#                               eta_timegrid_alpha_bamlss^2,
#                             index = ind,
#                             dX = NULL,
#                             deta = NULL,
#                             deta2 = NULL)



# Integration for multiple markers ----------------------------------------

nmarker <- 2

# Have the same function for the second longitudinal marker and the association
# This gives the function to integrate for marker-SPECIFIC effects:
# Intercept: \int exp(4+2x) dx -> 0.5*exp(4+2x)
# Time: \int x*exp(4+2x) dx -> 0.25*exp(2x+4)*(2x-1)
# => Hess: 0.5*exp(4+2x), 0.25*exp(2x+4)*(2x-1), 0.25*exp(2x+4)*(2x-1),
#          0.25*exp(2x+4)*(2x^2-2x+1)

# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(0.5*exp(4+2*t)-0.5*exp(4),
                     0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4)),
       hess_int = c(0.5*exp(4+2*t)-0.5*exp(4),
                    0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4),
                    0.25*exp(2*t+4)*(2*t-1) + 0.25*exp(4),
                    0.25*exp(2*t+4)*(2*t^2-2*t+1) - 0.25*exp(4)))
}

# Set up arguments for survint
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
int_new
int_func(s_times[1])


# This gives the function to integrate for marker-COMMON effects:
# Intercept: \int 2*exp(4+2x) dx -> exp(4+2x)
# Time: \int 2x*exp(4+2x) dx -> 0.5*exp(2x+4)*(2x-1)
# => Hess: exp(4+2x), 0.5*exp(2x+4)*(2x-1), 0.5*exp(2x+4)*(2x-1),
#          0.5*exp(2x+4)*(2x^2-2x+1)

# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(exp(4+2*t)-exp(4), 0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4)),
       hess_int = c(exp(4+2*t)-exp(4), 0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4),
                    0.5*exp(2*t+4)*(2*t-1) + 0.5*exp(4),
                    0.5*exp(2*t+4)*(2*t^2-2*t+1) - 0.5*exp(4)))
}

# Set up arguments for survint
x_Xgrid_re <- cbind(1, rep(int_times, nmarker))

int_new_re <- survint_gq(pred = "long",
                         pre_fac = exp(eta_gamma),
                         omega = exp(eta_timegrid),
                         int_fac = eta_timegrid_alpha,
                         int_vec = x_Xgrid_re,
                         weights = gq$weights,
                         survtime = s_times)
int_new_re
int_func(s_times[1])





