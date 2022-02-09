source("R/opt_updating.R")
source("R/survint.R")


# Set-Up for one longitudinal marker --------------------------------------

n <- 3
n_gq <- 7
subdivisions <- 250
nmarker <- 1

# Survival times
s_times <- seq(0, 1, length.out = (n+1))[-1]

# Integration GQ
gq <- statmod::gauss.quad(n = n_gq)
create_grid <- function(time) {
  time / 2 * gq$nodes + time / 2
}
int_times <- unlist(lapply(s_times, create_grid))

# Integration trapezoid
int_times_bamlss <- rep(1, n*subdivisions)
for (i in seq_len(n)) {
  int_times_bamlss[((i-1)*subdivisions + 1):(i*subdivisions)] <-
    seq(0, s_times[i], length.out = subdivisions)
}
width <- int_times_bamlss[(seq_len(n) - 1)*subdivisions + 2]




# Specify predictor functions ---------------------------------------------

lambda <- Vectorize(function (t) {
  -0.8*t
})
gamma <- function (x) {
  1*x
}
alpha <- Vectorize(function (x) {
  0.5*x
})
mu <- function (t) {
  1 + 1*t
}


# This gives the functions to integrate (lambda):
# exp(-0.8*t+0.5*0.3(1+t))t [SCORE]
# exp(-0.8*t+0.5*0.3(1+t))tt [HESS]
int_lambda_one <- function (t, x) {
  list(score_int = (10*exp((5*x+t*(-8+5*x))/10)*(-10+t*(-8+5*x)))/(8-5*x)^2 -
         (10*exp(5*x/10)*(-10))/(8-5*x)^2,
       hess_int = (10*exp((5*x+t*(-8+5*x))/10)*(200+t^2*(8-5*x)^2-20*t*
                                                  (-8+5*x)))/(-8+5*x)^3 -
         (10*exp(5*x/10)*(200))/(-8+5*x)^3)
}
int_lambda <- function (t, x) {
  do.call(Map, c(f = rbind, mapply(int_lambda_one, t = t, x = x, 
                                   SIMPLIFY = FALSE)))
}

# This gives the functions to integrate (mu - Intercept):
# exp(-0.8*t+0.5*x*(1+t))*(0.5*x) [SCORE]
# exp(-0.8*t+0.5*x*(1+t))*(0.5*x)^2 [HESS]
int_mu_one <- function (t, x) {
  list(score_int = c((5*exp(x/2+(t*(-8+5*x))/10)*x)/(-8+5*x) - 
                       (5*exp(x/2)*x)/(-8+5*x),
                     5*exp((5*x+t*(-8+5*x))/10)*x*(-10+t*(-8+5*x))/(8-5*x)^2 -
                       5*exp((5*x)/10)*x*(-10)/(8-5*x)^2),
       hess_int = c((5*exp(x/2+(t*(-8+5*x))/10)*x^2)/(2*(-8+5*x)) -
                      (5*exp(x/2)*x^2)/(2*(-8+5*x)),
                    (5*exp((5*x+t*(-8+5*x))/10)*x^2*(-10+t*(-8+5*x)))/
                      (2*(8-5*x)^2) - (5*exp(5*x/10)*x^2*(-10))/(2*(8-5*x)^2),
                    (5*exp((5*x+t*(-8+5*x))/10)*x^2*(-10+t*(-8+5*x)))/
                      (2*(8-5*x)^2) - (5*exp(5*x/10)*x^2*(-10))/(2*(8-5*x)^2),
                    (5*exp((5*x+t*(-8+5*x))/10)*x^2*
                       (200+t^2*(8-5*x)^2-20*t*(-8+5*x)))/(2*(-8+5*x)^3) - 
                      (5*exp(5*x/10)*x^2*200)/(2*(-8+5*x)^3)))
}
int_mu <- function (t, x) {
  do.call(Map, c(f = rbind, mapply(int_mu_one, t = t, x = x, 
                                   SIMPLIFY = FALSE)))
}


# Set-Up design matrices --------------------------------------------------

# GQ
eta_timegrid_lambda <- lambda(int_times)
eta_gamma <- gamma(c(2, 1, 0.5))
eta_timegrid_alpha <- alpha(rep(c(2, 1, 0.5), each = n_gq))
eta_timegrid_mu <- mu(int_times)
eta_timegrid_long <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu))
eta_timegrid <- eta_timegrid_lambda + eta_timegrid_long
Xgrid_lambda <- cbind(int_times)
X_gamma <- cbind(c(2, 1, 0.5))
Xgrid_mu <- cbind(1, int_times)

# Trapezoidal
eta_timegrid_lambda_bamlss <- matrix(lambda(int_times_bamlss), nrow = n,
                                     byrow = TRUE)
eta_timegrid_alpha_bamlss <- matrix(alpha(rep(c(2, 1, 0.5), 
                                              each = subdivisions)),
                                    nrow = n, byrow = TRUE)
eta_timegrid_mu_bamlss <- matrix(mu(int_times_bamlss), nrow = n, byrow = TRUE)
eta_timegrid_bamlss <- eta_timegrid_lambda_bamlss +
  eta_timegrid_alpha_bamlss * eta_timegrid_mu_bamlss
X_lambda <- cbind(int_times_bamlss)
X_mu <- cbind(1, int_times_bamlss)



# Lamba -------------------------------------------------------------------

# GQ
survint_gq(pred = "lambda", 
           pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid),
           int_vec = Xgrid_lambda, 
           weights = gq$weights,
           survtime = s_times)

# GQ using for-loop
(gq_lambda <- survint_gqFOR(pred = "lambda", 
                             pre_fac = exp(eta_gamma),
                             omega = exp(eta_timegrid),
                             int_vec = Xgrid_lambda, 
                             weights = gq$weights,
                             survtime = s_times))

# Trapezoidal
bamlss:::survint(X = X_lambda,
                 eta = exp(eta_timegrid_bamlss),
                 width = width,
                 gamma = exp(eta_gamma),
                 eta2 = NULL,
                 index = NULL,
                 dX = NULL,
                 deta = NULL,
                 deta2 = NULL)
list(score = colSums(gq_lambda$score_int), 
     hess = matrix(colSums(gq_lambda$hess_int), nrow = 1))

# True integral
int_lt <- int_lambda(s_times, eta_gamma)
list(score = exp(eta_gamma) %*% int_lt$score_int, 
     hess = matrix(exp(eta_gamma) %*% int_lt$hess_int, nrow = 1))


# Gamma -------------------------------------------------------------------

# GQ
survint_gq(pred = "gamma",
           pre_fac = exp(eta_gamma),
           pre_vec = X_gamma,
           omega = exp(eta_timegrid),
           weights = gq$weights,
           survtime = s_times)

# GQ using for-loop
gq_gamma <- survint_gqFOR(pred = "gamma",
                          pre_fac = exp(eta_gamma),
                          pre_vec = X_gamma,
                          omega = exp(eta_timegrid),
                          weights = gq$weights,
                          survtime = s_times)
colSums(gq_gamma$score_int)
colSums(gq_gamma$hess_int)

# Trapezoidal implementation in bamlss
int0 <- width * (0.5 * (exp(eta_timegrid_bamlss)[, 1] + 
                          exp(eta_timegrid_bamlss)[, subdivisions]) + 
                   apply(exp(eta_timegrid_bamlss)[, 2:(subdivisions - 1)], 1,
                         sum))
tint <- xhess0 <- 0
for(i in 1:n) {
  tint <- tint + exp(eta_gamma[i]) * X_gamma[i, , drop = FALSE] * int0[i]
  xhess0 <- xhess0 + exp(eta_gamma[i]) * X_gamma[i, ] %*% t(X_gamma[i, ]) *
    int0[i]
}
tint
xhess0



# Longitudinal ------------------------------------------------------------

# GQ
survint_gq(pred = "long", pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid),
           int_fac = eta_timegrid_alpha,
           int_vec = Xgrid_mu,
           weights = gq$weights,
           survtime = s_times)

# GQ using for-loop
(gq_mu <- survint_gqFOR(pred = "long", 
                        pre_fac = exp(eta_gamma),
                        omega = exp(eta_timegrid),
                        int_fac = eta_timegrid_alpha,
                        int_vec = Xgrid_mu, 
                        weights = gq$weights,
                        survtime = s_times))
list(score = colSums(gq_mu$score_int), 
     hess = matrix(colSums(gq_mu$hess_int), nrow = 2))

# Trapezoidal
bamlss:::survint(X = X_mu,
                 eta = exp(eta_timegrid_bamlss) * eta_timegrid_alpha_bamlss,
                 width = width,
                 gamma = exp(eta_gamma),
                 eta2 = exp(eta_timegrid_bamlss) * eta_timegrid_alpha_bamlss^2,
                 index = NULL,
                 dX = NULL,
                 deta = NULL,
                 deta2 = NULL)

# True integral
int_mt <- int_mu(s_times, eta_gamma)
exp(eta_gamma) * int_mt$score_int
list(score = exp(eta_gamma) %*% int_mt$score_int, 
     hess = matrix(exp(eta_gamma) %*% int_mt$hess_int, nrow = 2))




# Multiple longitudinal markers -------------------------------------------

nmarker <- 2
alpha_1 <- Vectorize(function (x) {
  0.5*x
})
alpha_2 <- Vectorize(function (x) {
  -3
})
mu_1 <- function (t) {
  1 + 1*t
}
mu_2 <- function (t) {
  0.5*t
}

# This gives the functions to integrate (mu1 - Intercept):
# exp(-0.8*t+0.5*x*(1+t)-3*0.5*t)*(0.5*x) [SCORE]
# exp(-0.8*t+0.5*x*(1+t)-3*0.5*t)*(0.5*x)^2 [HESS]
int_mu1_one <- function (t, x) {
  list(score_int = c((5*exp(x/2+(t*(-23+5*x))/10)*x)/(-23+5*x) - 
                       (5*exp(x/2)*x)/(-23+5*x),
                     5*exp((5*x+t*(-8+5*x))/10)*x*(-10+t*(-8+5*x))/(8-5*x)^2 -
                       5*exp((5*x)/10)*x*(-10)/(8-5*x)^2),
       hess_int = c((5*exp(x/2+(t*(-23+5*x))/10)*x^2)/(2*(-23+5*x)) -
                      (5*exp(x/2)*x^2)/(2*(-23+5*x)),
                    (5*exp((5*x+t*(-23+5*x))/10)*x^2*(-10+t*(-23+5*x)))/
                      (2*(23-5*x)^2) - (5*exp(5*x/10)*x^2*(-10))/(2*(23-5*x)^2),
                    (5*exp((5*x+t*(-23+5*x))/10)*x^2*(-10+t*(-23+5*x)))/
                      (2*(23-5*x)^2) - (5*exp(5*x/10)*x^2*(-10))/(2*(23-5*x)^2),
                    (5*exp((5*x+t*(-8+5*x))/10)*x^2*
                       (200+t^2*(8-5*x)^2-20*t*(-8+5*x)))/(2*(-8+5*x)^3) - 
                      (5*exp(5*x/10)*x^2*200)/(2*(-8+5*x)^3)))
}
int_mu1 <- function (t, x) {
  do.call(Map, c(f = rbind, mapply(int_mu1_one, t = t, x = x, 
                                   SIMPLIFY = FALSE)))
}


eta_timegrid_alpha_mul <- c(alpha_1(rep(c(2, 1, 0.5), each = n_gq)), 
                            alpha_2(rep(1, n*n_gq)))
eta_timegrid_mu_mul <- c(mu_1(int_times), mu_2(int_times))
eta_timegrid_long_mul <- rowSums(
  matrix(eta_timegrid_alpha_mul*eta_timegrid_mu_mul, ncol = nmarker, 
         nrow = n*n_gq))
eta_timegrid_mul <- eta_timegrid_lambda + eta_timegrid_long_mul
Xgrid_mu_mul <- rbind(Xgrid_mu, matrix(0, nrow = n*n_gq, ncol = 2))


# GQ
survint_gq(pred = "long", pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid_mul),
           int_fac = eta_timegrid_alpha_mul,
           int_vec = Xgrid_mu_mul,
           weights = gq$weights,
           survtime = s_times)
survint_gqFOR(pred = "long", pre_fac = exp(eta_gamma),
              omega = exp(eta_timegrid_mul),
              int_fac = eta_timegrid_alpha_mul,
              int_vec = Xgrid_mu_mul,
              weights = gq$weights,
              survtime = s_times)

# True integral
int_mu1t <- int_mu1(s_times, eta_gamma)
exp(eta_gamma) * int_mt$score_int
