load("inst/pd/221017_survint.Rdata")
source("R/survint.R")
source("R/compile.R")
compile_alex("workstation")

c_long <- survint_C(pred = pred, pre_fac = pre_fac, omega = omega, 
                    int_fac = int_fac, int_vec = int_vec, weights = weights,
                    survtime = survtime)
c_re <- survint_C(pred = "fpc_re", pre_fac = pre_fac, omega = omega, 
                  int_fac = int_fac, int_vec = int_vec, weights = weights,
                  survtime = survtime)
r_long <- survint_gq0(pred = pred, pre_fac = pre_fac, omega = omega, 
                      int_fac = int_fac, int_vec = int_vec, weights = weights,
                      survtime = survtime)
all.equal(colSums(r_long$score_int), c_long$score_int)
all.equal(colSums(r_long$hess_int), c_long$hess_int)
all.equal(c_long$score_int, c_re$score_int)
all.equal(c_long$hess_int, c(diag(c_re$hess_int)))

