load(file = "Fehlersuche/Unterschiede_survint/survint_bamlss.RData")
int <- bamlss:::survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                        eeta * eta_timegrid_alpha^2,
                        index = x$sparse.setup[["mu.matrix"]],
                        dX, NULL, NULL)
eta_bamlss <- eta
x_bamlss <- x

source("R/survint.R")
load(file = "Fehlersuche/Unterschiede_survint/survint_mjm.RData")
int_i <- survint_gq(pred = "long", pre_fac = exp(eta$gamma),
                    omega = exp(eta_timegrid),
                    int_fac = eta_timegrid_alpha, int_vec = x$Xgrid,
                    weights = gq_weights, survtime = survtime)
# Quite large computationally
# Or load result
# load("Fehlersuche/Unterschiede_survint/int_i.RData")


# Compare survival integral in gradient
int$grad
colSums(int_i$score_int)

# Compare survival integral in hesse
int$hess
matrix(colSums(int_i$hess_int), ncol = ncol(int_i$score_int))
