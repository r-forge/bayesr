
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

# Integration trapezoid
int_times_bamlss <- rep(1, n*subdivisions)
for (i in seq_len(n)) {
  int_times_bamlss[((i-1)*subdivisions + 1):(i*subdivisions)] <-
    seq(0, s_times[i], length.out = subdivisions)
}
width <- int_times_bamlss[(seq_len(n) - 1)*subdivisions + 2]

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
int_one <- function (t) {
  do.call(Map, c(f = rbind, lapply(t, int_func_one)))
}

# survint_gq -----------
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

(int_i_one <- survint_gq(pred = "long",
                         pre_fac = exp(eta_gamma),
                         omega = exp(eta_timegrid),
                         int_fac = eta_timegrid_alpha,
                         int_vec = x_Xgrid,
                         weights = gq$weights,
                         survtime = s_times))
int_gq_one <- list(score = colSums(int_i_one$score_int), 
                   hess = matrix(colSums(int_i_one$hess_int),nrow = 2))

# Compare to correct integrals
(int_t_one <- int_one(s_times))
int_true_one <- list(score = colSums(int_t_one$score_int), 
                     hess = matrix(colSums(int_t_one$hess_int),nrow = 2))


# bamlss ----------------
# Integration grid
eta_timegrid_lambda_bamlss <- matrix(lambda(rep(1, n*subdivisions)), nrow = n,
                                     byrow = TRUE)
eta_timegrid_alpha_bamlss <- matrix(alpha(rep(1, n*subdivisions)), nrow = n,
                                    byrow = TRUE)
eta_timegrid_mu_bamlss <- matrix(mu(int_times_bamlss), nrow = n, byrow = TRUE)
eta_timegrid_bamlss <- eta_timegrid_lambda_bamlss +
  eta_timegrid_alpha_bamlss * eta_timegrid_mu_bamlss
X <- cbind(1, int_times_bamlss)

(int_old_one <- bamlss:::survint(X = X,
                                 eta = exp(eta_timegrid_bamlss) *
                                   eta_timegrid_alpha_bamlss,
                                 width = width,
                                 gamma = exp(eta_gamma),
                                 eta2 = exp(eta_timegrid_bamlss) *
                                   eta_timegrid_alpha_bamlss^2,
                                 index = NULL,
                                 dX = NULL,
                                 deta = NULL,
                                 deta2 = NULL))

# Compare to correct integrals
int_true_one
int_gq_one



# Check the survint function for non-unit interval ------------------------

# Survival times
s_times_120 <- seq(0, 120, length.out = (n+1))[-1]
int_times_120 <- unlist(lapply(s_times_120, create_grid))
int_times_120_bamlss <- rep(1, n*subdivisions)
for (i in seq_len(n)) {
  int_times_120_bamlss[((i-1)*subdivisions + 1):(i*subdivisions)] <-
    seq(0, s_times_120[i], length.out = subdivisions)
}
width_120 <- int_times_120_bamlss[(seq_len(n) - 1)*subdivisions + 2]

# survint_gq -------
eta_timegrid_mu_120 <- mu(int_times_120)
eta_timegrid_long_120 <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha*eta_timegrid_mu_120))
eta_timegrid_120 <- eta_timegrid_lambda + eta_timegrid_long_120
x_Xgrid_120 <- cbind(1, int_times_120)

# bamlss -----------
eta_timegrid_mu_120_bamlss <- matrix(mu(int_times_120_bamlss), nrow = n, 
                                     byrow = TRUE)
eta_timegrid_120_bamlss <- eta_timegrid_lambda_bamlss +
  eta_timegrid_alpha_bamlss * eta_timegrid_mu_120_bamlss
X_120 <- cbind(1, int_times_120_bamlss)


(int_t_120 <- int_one(s_times_120))
(int_i_120 <- survint_gq(pred = "long",
                         pre_fac = exp(eta_gamma),
                         omega = exp(eta_timegrid_120),
                         int_fac = eta_timegrid_alpha,
                         int_vec = x_Xgrid_120,
                         weights = gq$weights,
                         survtime = s_times_120))
(int_old_120 <- bamlss:::survint(X = X_120,
                                 eta = exp(eta_timegrid_120_bamlss) *
                                   eta_timegrid_alpha_bamlss,
                                 width = width_120,
                                 gamma = exp(eta_gamma),
                                 eta2 = exp(eta_timegrid_120_bamlss) *
                                   eta_timegrid_alpha_bamlss^2,
                                 index = NULL,
                                 dX = NULL,
                                 deta = NULL,
                                 deta2 = NULL))
(int_gq_120 <- list(score = colSums(int_i_120$score_int), 
                    hess = matrix(colSums(int_i_120$hess_int),nrow = 2)))
(int_true_120 <- list(score = colSums(int_t_120$score_int), 
                      hess = matrix(colSums(int_t_120$hess_int),nrow = 2)))

# ------LARGE DIFFERENCES--------
# Both underestimate the integral but GQ is closer to the truth
abs(int_true_120$score - int_gq_120$score) - 
  abs(int_true_120$score - int_old_120$grad) < 0



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
int_one_dec <- function (t) {
  do.call(Map, c(f = rbind, lapply(t, int_func_one_dec)))
}



# Set up arguments for survint --------------------------------------
# survint_gq ------
eta_timegrid_lambda_dec <- lambda_dec(int_times)
eta_gamma_dec <- gamma_dec(rep(1, n))
eta_timegrid_alpha_dec <- alpha_dec(rep(1, n*n_gq))
eta_timegrid_mu_dec <- mu_dec(int_times)
eta_timegrid_long_dec <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha_dec*eta_timegrid_mu_dec))
eta_timegrid_dec <- eta_timegrid_lambda_dec + eta_timegrid_long_dec

# bamlss ----------
eta_timegrid_lambda_decbam <- matrix(lambda_dec(rep(1, n*subdivisions)),
                                     nrow = n, byrow = TRUE)
eta_timegrid_alpha_decbam <- matrix(alpha_dec(rep(1, n*subdivisions)), nrow = n,
                                    byrow = TRUE)
eta_timegrid_mu_decbam <- matrix(mu_dec(int_times_bamlss), nrow = n,
                                 byrow = TRUE)
eta_timegrid_decbam <- eta_timegrid_lambda_decbam +
  eta_timegrid_alpha_decbam * eta_timegrid_mu_decbam


# Compare the different implementations
int_i_dec <- survint_gq(pred = "long",
                        pre_fac = exp(eta_gamma_dec),
                        omega = exp(eta_timegrid_dec),
                        int_fac = eta_timegrid_alpha_dec,
                        int_vec = x_Xgrid,
                        weights = gq$weights,
                        survtime = s_times)
int_t_dec <- int_one_dec(s_times)
(int_gq_dec <- list(score = colSums(int_i_dec$score_int), 
                    hess = matrix(colSums(int_i_dec$hess_int),nrow = 2)))
(int_true_dec <- list(score = colSums(int_t_dec$score_int), 
                      hess = matrix(colSums(int_t_dec$hess_int),nrow = 2)))
(int_old_one <- bamlss:::survint(X = X,
                                 eta = exp(eta_timegrid_decbam) *
                                   eta_timegrid_alpha_decbam,
                                 width = width,
                                 gamma = exp(eta_gamma_dec),
                                 eta2 = exp(eta_timegrid_decbam) *
                                   eta_timegrid_alpha_decbam^2,
                                 index = NULL,
                                 dX = NULL,
                                 deta = NULL,
                                 deta2 = NULL))





# Now with non-unit interval ----------------------------------------------


# survint_gq -------
eta_timegrid_mu_dec_120 <- mu_dec(int_times_120)
eta_timegrid_long_dec_120 <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha_dec*eta_timegrid_mu_dec_120))
eta_timegrid_dec_120 <- eta_timegrid_lambda_dec + eta_timegrid_long_dec_120

# bamlss -----------
eta_timegrid_mu_120_decbam <- matrix(mu_dec(int_times_120_bamlss), nrow = n, 
                                     byrow = TRUE)
eta_timegrid_120_decbam <- eta_timegrid_lambda_decbam +
  eta_timegrid_alpha_decbam * eta_timegrid_mu_120_decbam


(int_t_dec_120 <- int_one_dec(s_times_120))
(int_i_dec_120 <- survint_gq(pred = "long",
                             pre_fac = exp(eta_gamma_dec),
                             omega = exp(eta_timegrid_dec_120),
                             int_fac = eta_timegrid_alpha_dec,
                             int_vec = x_Xgrid_120,
                             weights = gq$weights,
                             survtime = s_times_120))
(int_gq_dec_120 <- list(score = colSums(int_i_dec_120$score_int), 
                        hess = matrix(colSums(int_i_dec_120$hess_int),
                                      nrow = 2)))
(int_true_dec <- list(score = colSums(int_t_dec_120$score_int), 
                      hess = matrix(colSums(int_t_dec_120$hess_int),nrow = 2)))
(int_old_120 <- bamlss:::survint(X = X_120,
                                 eta = exp(eta_timegrid_120_decbam) *
                                   eta_timegrid_alpha_decbam,
                                 width = width_120,
                                 gamma = exp(eta_gamma_dec),
                                 eta2 = exp(eta_timegrid_120_decbam) *
                                   eta_timegrid_alpha_decbam^2,
                                 index = NULL,
                                 dX = NULL,
                                 deta = NULL,
                                 deta2 = NULL))



# Complex functions with covariates ---------------------------------------

# constant baseline hazard, constant association
lambda_comp <- Vectorize(function (t) {
  log(2*t+1)
})
gamma_comp <- function (x) {
  1*x
}
alpha_comp <- Vectorize(function (x) {
  0.5*x
})
mu_comp <- function (t) {
  1 + 1*t
}

# This gives the functions to integrate:
# Intercept: exp(log(2*x+1)+0.5*b*(1+x))*(0.5*b) [SCORE]
# exp(log(2*x+1)+0.5*b*(1+x))*(0.5*b)^2 [HESS]
int_func_comp <- function (x, b) {
  list(score_int = c(exp(0.5*b*(x+1))*(2*b*x+b-4) / b - exp(0.5*b)*(b-4) / b,
                     exp(0.5*b*(x+1))*(b^2*x*(2*x+1)-2*b*(4*x+1)+16)/b^2 -
                       exp(0.5*b)*(-2*b+16)/b^2),
       hess_int = c(0.5*exp(0.5*b*(x+1))*(2*b*x+b-4) -
                      0.5*exp(0.5*b)*(b-4),
                    (2*b^2*x^2+(b^2-8*b)*x-2*b+16)*exp(0.5*b*(x+1))/(2*b) -
                      (-2*b+16)*exp(0.5*b)/(2*b),
                    (2*b^2*x^2+(b^2-8*b)*x-2*b+16)*exp(0.5*b*(x+1))/(2*b) -
                      (-2*b+16)*exp(0.5*b)/(2*b),
                    ((2*b^3*x^3+(b^3-12*b^2)*x^2+(48*b-4*b^2)*x+8*b-96)*
                        exp((b*x)/2+b/2))/(2*b^2) - 
                      ((8*b-96)*exp(b/2))/(2*b^2)))
}
int_comp <- function (t, x) {
  do.call(Map, c(f = rbind, mapply(int_func_comp, x = t, b = x, 
                                   SIMPLIFY = FALSE)))
}


# survint_gq ------------
eta_timegrid_lambda_comp <- lambda_comp(int_times)
eta_gamma_comp <- gamma(c(2, 1, 0.5))
eta_timegrid_alpha_comp <- alpha_comp(rep(c(2, 1, 0.5), each = n_gq))
eta_timegrid_mu_comp <- mu_comp(int_times)
eta_timegrid_long_comp <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda_comp)) %*%
    (eta_timegrid_alpha_comp*eta_timegrid_mu_comp))
eta_timegrid_comp <- eta_timegrid_lambda_comp + eta_timegrid_long_comp

# bamlss ----------------
eta_timegrid_lambda_bamcomp <- matrix(lambda_comp(int_times_bamlss), 
                                      nrow = n, byrow = TRUE)
eta_timegrid_alpha_bamcomp <- matrix(alpha_comp(rep(c(2, 1, 0.5), 
                                                    each = subdivisions)),
                                     nrow = n, byrow = TRUE)
eta_timegrid_mu_bamcomp <- matrix(mu_comp(int_times_bamlss), nrow = n,
                                  byrow = TRUE)
eta_timegrid_bamcomp <- eta_timegrid_lambda_bamcomp +
  eta_timegrid_alpha_bamcomp * eta_timegrid_mu_bamcomp


# Compare the different implementations
(int_i_comp <- survint_gq(pred = "long",
                         pre_fac = exp(eta_gamma_comp),
                         omega = exp(eta_timegrid_comp),
                         int_fac = eta_timegrid_alpha_comp,
                         int_vec = x_Xgrid,
                         weights = gq$weights,
                         survtime = s_times))
int_gq_comp <- list(score = colSums(int_i_comp$score_int), 
                   hess = matrix(colSums(int_i_comp$hess_int),nrow = 2))
(int_t_comp <- int_comp(s_times, eta_gamma_comp))
int_true_comp <- list(score = exp(eta_gamma_comp) %*% int_t_comp$score_int, 
                     hess = matrix(exp(eta_gamma_comp) %*% int_t_comp$hess_int,
                                   nrow = 2))
(int_old_comp <- bamlss:::survint(X = X,
                                  eta = exp(eta_timegrid_bamcomp) *
                                    eta_timegrid_alpha_bamcomp,
                                  width = width,
                                  gamma = exp(eta_gamma_comp),
                                  eta2 = exp(eta_timegrid_bamcomp) *
                                    eta_timegrid_alpha_bamcomp^2,
                                  index = NULL,
                                  dX = NULL,
                                  deta = NULL,
                                  deta2 = NULL))
int_true_comp
int_gq_comp




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

survint_gq(pred = "long",
           pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid),
           int_fac = eta_timegrid_alpha,
           int_vec = x_Xgrid,
           weights = gq$weights,
           survtime = s_times)
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

survint_gq(pred = "long",
           pre_fac = exp(eta_gamma),
           omega = exp(eta_timegrid),
           int_fac = eta_timegrid_alpha,
           int_vec = x_Xgrid_re,
           weights = gq$weights,
           survtime = s_times)
int_func(s_times[1])



# Integration multiple markers, decreasing hazard, non-unit -----------------


# Have the same function for the second longitudinal marker and the association
# This gives the function to integrate for marker-SPECIFIC effects:
# Intercept: \int exp(-4-0.02x) dx -> -50*exp(-4-0.02x)
# Time: \int x*exp(-4-0.02x) dx -> 50*exp(-0.02x-4)*(x+50)
# => Hess: 0.5*exp(4+2x), 0.25*exp(2x+4)*(2x-1), 0.25*exp(2x+4)*(2x-1),
#          0.25*exp(2x+4)*(2x^2-2x+1)

# Correct integrals would be
int_func <- function (t) {
  list(score_int = c(-50*exp(-4-0.02*t)+50*exp(-4),
                     -50*exp(-0.02*t-4)*(t+50) + 2500*exp(-4)),
       hess_int = c(-50*exp(-4-0.02*t)+50*exp(-4),
                    -50*exp(-0.02*t-4)*(t+50) + 2500*exp(-4),
                    -50*exp(-0.02*t-4)*(t+50) + 2500*exp(-4),
                    -50*exp(-0.02*t-4)*(t^2+100*t+5000) + 250000*exp(-4)))
}


# Set up arguments for survint
eta_timegrid_lambda_dec <- lambda_dec(int_times)
eta_gamma_dec <- gamma_dec(rep(1, n))
eta_timegrid_alpha_dec <- alpha_dec(rep(rep(1, n*n_gq), nmarker))
eta_timegrid_mu_dec <- mu_dec(rep(int_times_120, nmarker))
eta_timegrid_long_dec <- drop(
  t(rep(1, nmarker)) %x% diag(length(eta_timegrid_lambda)) %*%
    (eta_timegrid_alpha_dec*eta_timegrid_mu_dec))
eta_timegrid_dec <- eta_timegrid_lambda_dec + eta_timegrid_long_dec
x_Xgrid_120 <- rbind(cbind(1, int_times_120), 
                     cbind(0, rep(0, length(int_times))))

survint_gq(pred = "long",
           pre_fac = exp(eta_gamma_dec),
           omega = exp(eta_timegrid_dec),
           int_fac = eta_timegrid_alpha_dec,
           int_vec = x_Xgrid_120,
           weights = gq$weights,
           survtime = s_times_120)
int_func(s_times_120[1])

