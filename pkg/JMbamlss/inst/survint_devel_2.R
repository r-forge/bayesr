source("R/opt_updating.R")
source("R/survint.R")
source("R/compile.R")
compile_alex()

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


# Set-Up design matrices --------------------------------------------------


# Example
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


# Survint Calls -----------------------------------------------------------

survint_gq0(pred = "lambda", 
            pre_fac = exp(eta_gamma),
            omega = exp(eta_timegrid),
            int_vec = Xgrid_lambda, 
            weights = gq$weights,
            survtime = s_times)

survint_C(pred = "gamma",
          pre_fac = exp(eta_gamma),
          pre_vec = X_gamma,
          omega = exp(eta_timegrid),
          weights = gq$weights,
          survtime = s_times)

survint_gq0(pred = "gamma",
             pre_fac = exp(eta_gamma),
             pre_vec = X_gamma,
             omega = exp(eta_timegrid),
             weights = gq$weights,
             survtime = s_times)
survint_alex(pred = "gamma",
             pre_fac = exp(eta_gamma),
             pre_vec = X_gamma,
             omega = exp(eta_timegrid),
             weights = gq$weights,
             survtime = s_times)

survint_alex(pred = "gamma",
             pre_fac = exp(eta_gamma),
             pre_vec = X_gamma,
             omega = exp(eta_timegrid),
             weights = gq$weights,
             survtime = s_times)

